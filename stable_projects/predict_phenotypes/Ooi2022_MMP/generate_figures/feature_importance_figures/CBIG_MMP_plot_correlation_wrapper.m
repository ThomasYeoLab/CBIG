function CBIG_MMP_plot_correlation_wrapper
% CBIG_MMP_plot_correlation_wrapper
% This wrapper function calculates the correlation between feature
% importance values for each single-type-feature model for each modality.
% It additionally calculates the correlation between FC models for the
% cognition component between the HCP and ABCD.
%
% Prerequisites: Please run CBIG_MMP_<dataset>_plot_importance_wrapper.m
% first to get the collated feature importance mat files first.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set path to results
HCP_int_dir = fullfile('/home', 'leon_ooi', 'storage','Multimodal_prediction_project', ...
    'replication', 'HCP', 'output', 'interpretation');
ABCD_int_dir = fullfile('/home', 'leon_ooi', 'storage','Multimodal_prediction_project', ...
    'replication', 'ABCD', 'output', 'interpretation');

% set variables for behaviors and models
modalities = {'t1' 'tbss' 'tractography' 'fmri'};

% set HCP variables
HCP_outstems.t1 = {'Thickness' 'Area' 'Volume'};
HCP_outstems.tbss = {'FA' 'MD' 'AD' 'RD' 'OD' 'ICVF' 'ISOVF'};
HCP_outstems.tractography = {'FA' 'MD' 'AD' 'RD' 'OD' 'ICVF' 'ISOVF' ...
    'Stream count (log)' 'Stream length' };
HCP_outstems.fmri = {'Resting' 'Social' 'Gambling' ...
    'Language' 'Working Memory' 'Motor'};
HCP_cog = 60;
HCP_dir = fullfile(HCP_int_dir, 'mean_results');

% set ABCD variables
ABCD_outstems.t1 = {'Thickness' 'Area' 'Volume'};
ABCD_outstems.tbss = {'FA' 'MD' 'AD' 'RD' 'OD' 'ICVF' 'ISOVF'};
ABCD_outstems.tractography = {'FA' 'MD' 'AD' 'RD' 'OD' 'ICVF' 'ISOVF' ...
    'Stream count (log)' 'Stream length' };
ABCD_outstems.fmri = {'Resting' 'N-back' 'MID' 'SST'};

ABCD_cog = 37;
ABCD_dir = fullfile(ABCD_int_dir, 'mean_results');

% start plotting correlation
for n = 1:length(modalities)
    % get current modality
    curr_mod = modalities{n};
    HCP_names = HCP_outstems.(curr_mod);
    ABCD_names = ABCD_outstems.(curr_mod);
    % load reordered modality mean
    if strcmp(curr_mod, 'tbss')
        % read in tbss individual files for tract-averaged values
        for m = 1:length(HCP_outstems.tbss)
            curr_feat = HCP_outstems.tbss{m};
            % cognition
            HCP_tbss = load(fullfile(HCP_int_dir, strcat('features_tbss_', curr_feat), 'cog_fi_unsorted.mat'));
            ABCD_tbss = load(fullfile(ABCD_int_dir, strcat('features_tbss_', curr_feat), 'cog_fi_unsorted.mat'));
            HCP_imp.mod_mean_importance(m,:,60) = HCP_tbss.all_loadings(:,2);
            ABCD_imp.mod_mean_importance(m,:,37) = ABCD_tbss.all_loadings(:,2);
            % dissatisfaction / personality
            HCP_tbss = load(fullfile(HCP_int_dir, strcat('features_tbss_', curr_feat), 'dis_fi_unsorted.mat'));
            ABCD_tbss = load(fullfile(ABCD_int_dir, strcat('features_tbss_', curr_feat), 'pers_fi_unsorted.mat'));
            HCP_imp.mod_mean_importance(m,:,59) = HCP_tbss.all_loadings(:,2);
            ABCD_imp.mod_mean_importance(m,:,39) = ABCD_tbss.all_loadings(:,2);
            % emotion / mh
            HCP_tbss = load(fullfile(HCP_int_dir, strcat('features_tbss_', curr_feat), 'emo_fi_unsorted.mat'));
            ABCD_tbss = load(fullfile(ABCD_int_dir, strcat('features_tbss_', curr_feat), 'mh_fi_unsorted.mat'));
            HCP_imp.mod_mean_importance(m,:,61) = HCP_tbss.all_loadings(:,2);
            ABCD_imp.mod_mean_importance(m,:,38) = ABCD_tbss.all_loadings(:,2);
        end
    else
        ABCD_imp = load(fullfile(ABCD_dir,curr_mod,strcat(curr_mod, '_mean_imp.mat')));
        HCP_imp = load(fullfile(HCP_dir,curr_mod,strcat(curr_mod, '_mean_imp.mat')));
    end
    
    % reorder HCP fMRI so that it matches tasks in ABCD
    if strcmp(curr_mod, 'fmri')
        HCP_reorder = [1 5 3 6 2 4];
        HCP_imp.mod_mean_importance = HCP_imp.mod_mean_importance(HCP_reorder, :,:);
    end
    
    % Plot correlation across all behaviours
    plot_all_corrmat(HCP_imp.mod_mean_importance, ABCD_imp.mod_mean_importance, ...
        fullfile(HCP_dir,curr_mod), 'fullcorr')
    
    % Plot correlation across datasets for cognition
    HCP_cog_vec = HCP_imp.mod_mean_importance(:,:,HCP_cog)';
    ABCD_cog_vec = ABCD_imp.mod_mean_importance(:,:,ABCD_cog)';
    if size(HCP_cog_vec,1) ~= size(ABCD_cog_vec,1)
        fprintf('Features between HCP and ABCD have different lengths, skipping \n')
        continue
    end
    plot_within_between_corrmat(HCP_cog_vec, ABCD_cog_vec, ...
        fullfile(HCP_dir,curr_mod), "Cognition")
end
end

function plot_all_corrmat(HCP_behav, ABCD_behav, outdir, prefix)
% this function plots correlation of feature importance for each behaviour
% for predictive models of the same modality for the HCP and ABCD
% separately.
%
% Input:
% - HCP_behav
% A #models x #features x #behaviours mat file with importance values from
% each model from a specified modality for the HCP.
%
% - ABCD_behav
% A #models x #features x #behaviours mat file with importance values from
% each model from a specified modality for the ABCD.
%
% - outdir
% Directory for results to be saved.
%
% - prefix
% Prefix for the plots.

% plot for HCP: in the order of cognition, dissatisfaction, emotion
HCP_idx = [60 59 61];
% arrange plot for each behavior
plotsz = size(HCP_behav,1);
for n = 1:3
    idx_1 = HCP_idx(n);
    for m = 1:3
        idx_2 = HCP_idx(m);
        if n == m % plot values between each behavior
            pbox_1 = (n-1)*plotsz + 1 : n*plotsz; % allocate location in plot
            HCP_corr(pbox_1, pbox_1) = CBIG_corr(HCP_behav(:,:,idx_1)', HCP_behav(:,:,idx_1)');
        else % plot values across behaviors
            pbox_1 = (n-1)*plotsz + 1 : n*plotsz; % allocate location in plot
            pbox_2 = (m-1)*plotsz + 1 : m*plotsz;
            HCP_corr(pbox_2, pbox_1) = CBIG_corr(HCP_behav(:,:,idx_1)', HCP_behav(:,:,idx_2)');
            HCP_corr(pbox_1, pbox_2) = CBIG_corr(HCP_behav(:,:,idx_1)', HCP_behav(:,:,idx_2)')';
        end
    end
end

% generate image
imagesc(HCP_corr);
colormap(CBIG_MMP_GenerateColorscale('blue_brown', outdir))
% generate lines
[xline, yline, ymaj] = generateline(size(HCP_corr,1));
% lines for minor grids
minor_grid = 1:size(HCP_corr,1)-1;
patch(xline(:,minor_grid), yline(:,minor_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
% lines for major grids (separate behaviors)
major_grid = size(HCP_behav,1):size(HCP_behav,1):size(HCP_corr,1)-1;
patch(xline(:,major_grid), ymaj(:,major_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 1.5,'EdgeAlpha', 0.9);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 1.5,'EdgeAlpha', 0.9);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
% save figure and values
saveas(gcf, fullfile(outdir, strcat('HCP_', prefix, '.png'))); close(gcf);
csvwrite(fullfile(outdir, strcat('HCP_', prefix, '.csv')), HCP_corr)

% plot for HCP: in the order of cognition, personality and mental health
ABCD_idx = [37 39 38];        idx_1 = ABCD_idx(n);
% arrange plot for each behavior
plotsz = size(ABCD_behav,1);
for n = 1:3
    idx_1 = ABCD_idx(n);
    for m = 1:3
        idx_2 = ABCD_idx(m);
        if n == m % plot values between each behavior
            pbox_1 = (n-1)*plotsz + 1 : n*plotsz; % allocate location in plot
            ABCD_corr(pbox_1, pbox_1) = CBIG_corr(ABCD_behav(:,:,idx_1)', ABCD_behav(:,:,idx_1)');
        else % plot values across behaviors
            pbox_1 = (n-1)*plotsz + 1 : n*plotsz; % allocate location in plot
            pbox_2 = (m-1)*plotsz + 1 : m*plotsz;
            ABCD_corr(pbox_2, pbox_1) = CBIG_corr(ABCD_behav(:,:,idx_1)', ABCD_behav(:,:,idx_2)');
            ABCD_corr(pbox_1, pbox_2) = CBIG_corr(ABCD_behav(:,:,idx_1)', ABCD_behav(:,:,idx_2)')';
        end
    end
end

% generate image
imagesc(ABCD_corr);
colormap(CBIG_MMP_GenerateColorscale('blue_brown', outdir))
% generate lines
[xline, yline, ymaj] = generateline(size(ABCD_corr,1));
% lines for minor grids
minor_grid = 1:size(ABCD_corr,1)-1;
patch(xline(:,minor_grid), yline(:,minor_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
% lines for major grids (separate behaviors)
major_grid = size(ABCD_behav,1):size(ABCD_behav,1):size(ABCD_corr,1)-1;
patch(xline(:,major_grid), ymaj(:,major_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 1.5,'EdgeAlpha', 0.9);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 1.5,'EdgeAlpha', 0.9);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
% save figure and values
saveas(gcf, fullfile(outdir, strcat('ABCD_', prefix, '.png'))); close(gcf);
csvwrite(fullfile(outdir, strcat('ABCD_', prefix, '.csv')), ABCD_corr);

end

function plot_within_between_corrmat(HCP_behav, ABCD_behav, ...
    outdir, prefix)
% this function plots correlation of feature importance for each
% behaviour for predictive models of the same modality between the HCP
% and ABCD
%
% Input:
% - HCP_behav
% A #models x #features x #behaviours mat file with importance values from
% each model from a specified modality for the HCP.
%
% - ABCD_behav
% A #models x #features x #behaviours mat file with importance values from
% each model from a specified modality for the ABCD.
%
% - outdir
% Directory for results to be saved.
%
% - prefix
% Prefix for the plots.

% plot correlation matrix for HCP-ABCD
HCP_ABCD_corr = CBIG_corr(HCP_behav, ABCD_behav);
imagesc(HCP_ABCD_corr)

% generate image
colormap(CBIG_MMP_GenerateColorscale('blue_brown', outdir))
caxis([-0.8 0.8])
% generate lines
[xline, yline, ymaj] = generateline(size(HCP_ABCD_corr,1));
% lines for minor grids
minor_grid = 1:size(HCP_ABCD_corr,1)-1;
patch(xline(:,minor_grid), yline(:,minor_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', ...
    'edgecolor', 'k', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
% save figure
saveas(gcf, fullfile(outdir, strcat('HCP_ABCD_', prefix, '.png'))); close(gcf);
end

function [x, y, ymaj] = generateline(n)
% this function generates lines for plots of given size n

% define x and y lines
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
% For short grid (same as the figure boundary) ymaj = [ 0.5 n+0.5 nan]';
ymaj = repmat(ymaj,1,(n-1));

end
