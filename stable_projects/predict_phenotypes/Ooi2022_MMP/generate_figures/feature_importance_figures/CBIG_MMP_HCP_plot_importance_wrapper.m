function CBIG_MMP_HCP_plot_importance_wrapper(int_dir)
% CBIG_MMP_HCP_plot_importance_wrapper
% This is a wrapper function that plots the feature importance for all 
% all KRR models for the HCP.
%
% Please specify the directory where your interpretation results are saved.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Add path to result: Change for your own directory (EXAMPLE BELOW)
%int_dir = fullfile('/home', 'leon_ooi', 'storage', 'Multimodal_prediction_project', ...
%    'replication', 'HCP', 'output', 'interpretation');

% Initialize modalities and feature names features
modalities = {'t1' 'tbss' 'tractography' 'fmri'};

all_outstems.t1 = {'features_ct' 'features_ca' 'features_cv'};
all_outstems.tbss = {'features_tbss_FA' 'features_tbss_MD' ...
    'features_tbss_AD' 'features_tbss_RD' 'features_tbss_OD' ...
    'features_tbss_ICVF' 'features_tbss_ISOVF'};
all_outstems.tractography = {'features_schaefer_FA' 'features_schaefer_MD' ...
    'features_schaefer_AD' 'features_schaefer_RD' 'features_schaefer_OD' ...
    'features_schaefer_ICVF' 'features_schaefer_ISOVF' ...
    'features_schaefer_streamcount_log' 'features_schaefer_streamlen' };
all_outstems.fmri = {'features_rs' 'features_social' 'features_gamb' ...
    'features_lang' 'features_wm' 'features_motor'};

% plot function for each modality
for n = 1:length(modalities)
    curr_mod = all_outstems.(modalities{n});
    clear mod_mean_importance
    for m = 1:length(curr_mod)
        % plot for each feature
        curr_outstem = curr_mod{m};
        % pass flag to rearrange ordering for networks generated in 17 and
        % 7 network ordering
        if strcmp(modalities{n}, 'tractography')
            network_rearr = 17;
        elseif strcmp(modalities{n}, 'fmri')
            if strcmp(curr_outstem, 'features_rs')
                network_rearr = 17;
            else
                network_rearr = 7;
            end
        else
            network_rearr = 0;
        end
        load(fullfile(int_dir, curr_outstem , 'cov_mat_mean.mat'))
        mod_mean_importance(m,:,:) = modality_plot(modalities{n}, ...
            cov_mat_mean, int_dir, curr_outstem, network_rearr); 
    end
    
    % save collated results
    mean_dir = fullfile(int_dir, 'mean_results');
    curr_outstem = modalities{n};
    if ~exist(fullfile(mean_dir,curr_outstem))
        mkdir(fullfile(mean_dir,curr_outstem))
    end
    % save reordered modality mean
    save(fullfile(mean_dir,curr_outstem,strcat(curr_outstem, '_mean_imp.mat')), ....
        'mod_mean_importance');
end
end

function mean_importance = modality_plot(modality, cov_mat_mean, int_dir, ...
    curr_outstem, network_rearr)
% this function defines score indices to be plotted for feature importance
% and plots it in the appropriate way for each modality. 
%
% Input:
% - modality
% A string specifying which modality to plot. Can be "t1", "tbss", 
% "tractography" or "fmri".
%
% - cov_mat_mean
% A #folds x #features x #behaviours matrix. Ordering of behaviours should 
% correspond to ordering specified in behaviour file provided for the paper.
%
% - int_dir
% Interpretation directory where feature importance results are saved.
% Images will be saved there as well.
%
% - curr_outstem
% Outstem for current model feature importance is being plotted for.
%
% - network_rearr
% Specify ordering in which current network features are arranged in. Only
% works for tractography and fmri results. "7" will convert Schaefer_Yeo7
% ordering to Schaefer_Kong2022 ordering. "17" converts Schaefe_Yeo17
% ordering to Schaefer_Kong2022 ordering. "0" indicates no reordering to be
% done.
%
% Output: 
% - mean_importance
% Mean importance values with updated ordering based on "network_arr".

% define score indices: Ensure that this ordering is correct for you
cog_idx = 60;
dis_idx = 59;
emo_idx = 61;    
prefixes = {'cog' 'dis' 'emo'};

% calculate mean importance
mean_importance = squeeze(mean(cov_mat_mean,1));
imp_cog = mean_importance(:,cog_idx);
imp_dis = mean_importance(:,dis_idx);
imp_emo = mean_importance(:,emo_idx);
        
% choose which modality to plot
switch modality
    case "t1"
        fprintf("Plot T1\n")
        % cognition
        prefix = prefixes{1};
        CBIG_MMP_plot_cortical(imp_cog, fullfile(int_dir, curr_outstem), ...
            prefix)
        % dissatisfaction
        prefix = prefixes{2};
        CBIG_MMP_plot_cortical(imp_dis, fullfile(int_dir, curr_outstem), ...
            prefix)
        % emotion
        prefix = prefixes{3};
        CBIG_MMP_plot_cortical(imp_emo, fullfile(int_dir, curr_outstem), ...
            prefix)
    case "tbss"
        fprintf("Plot TBSS\n")
        % additionally need to define path to skeleton mask
        skeleton_mask = fullfile('/home', 'leon_ooi', 'storage', 'Multimodal_prediction_project', ...
            'replication', 'HCP', 'input', 'HCP_processed_data', 'mean_FA_skeleton_mask.nii.gz');
        % cognition
        prefix = prefixes{1};
        CBIG_MMP_plot_tbss(imp_cog, fullfile(int_dir, curr_outstem), prefix, ...
            skeleton_mask)
        % dissatisfaction
        prefix = prefixes{2};
        CBIG_MMP_plot_tbss(imp_dis, fullfile(int_dir, curr_outstem), prefix, ...
            skeleton_mask)
        % emotion
        prefix = prefixes{3};
        CBIG_MMP_plot_tbss(imp_emo, fullfile(int_dir, curr_outstem), prefix, ...
            skeleton_mask)
    case {"tractography", "fmri"}
        fprintf("Plot ROI2ROI networks\n")
        % reorder importances to Schaefer_Kong2022 ordering
        if network_rearr ~= 0
            fprintf("\t reordering networks...\n")
            for behav = 1:size(mean_importance, 2)
                mean_importance(:,behav) = CBIG_MMP_reorder_imp(...
                    mean_importance(:,behav), network_rearr);
            end
            imp_cog = mean_importance(:,cog_idx);
            imp_dis = mean_importance(:,dis_idx);
            imp_emo = mean_importance(:,emo_idx);
        end
        % cognition
        prefix = prefixes{1};
        CBIG_MMP_plot_ROI2ROI(imp_cog, fullfile(int_dir, curr_outstem), ...
            prefix);
        % dissatisfaction
        prefix = prefixes{2};
        CBIG_MMP_plot_ROI2ROI(imp_dis, fullfile(int_dir, curr_outstem), ...
            prefix);
        % emotion
        prefix = prefixes{3};
        CBIG_MMP_plot_ROI2ROI(imp_emo, fullfile(int_dir, curr_outstem), ...
            prefix);
end
end