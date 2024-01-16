function acc_vec = CBIG_HCP_read_model_results_example(outstem, outdir, N_seeds, N_folds, behav_ind, metric, store_new)

% acc_vec = CBIG_MMP_HCP_read_model_results_example(outstem, outdir, N_seeds, N_folds, behav_ind, metric, store_new)
%
% This function reads the results of KRR, LRR, Elasticnet, multiKRR and stacking. The result is
% saved in a matrix of #behav x #splits. Results are saved for LRR, Elasticnet and stacking since
% computational time for these are high.
%
% Inputs:
%   - outstem
%     Outstem of the file to be read. Should be in the format of
%     `regression_feature`, derived from any of the regression models
%     (e.g. CBIG_MMP_HCP_KRR_example.m).
%
%   - outdir
%     Full path of directory where the results from the regression models
%     are saved.
%
%   - N_seeds
%     Number of splits in the regression model.
%
%   - N_folds
%     Number of folds in the regression model.
%
%   - behav_ind
%     A vector corresponding to the behavioural indices to be extracted.
%
%   - metric
%     Metric to be read. Can be chosen from
%     {'corr', 'COD', 'predictive_COD', 'MAE' 'MAE_norm', 'MSE', 'MSE_norm'}.
%
%   - store_new
%     A logical. Flag to save a new mat file of acc_vec.
%
% Outputs:
%   - acc_vec
%     A matrix of accuracies in the dimensions of #behav x #splits.
%
% Written by Leon_Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

split_outstem = strsplit(outstem,'_');
model = split_outstem{1};
N_behav = length(behav_ind);

%% read results based on model type
switch model
    %% first level models
    case 'KRR'
        for seed_num = 1:N_seeds
            seed_name = strcat('seed_', num2str(seed_num));
            tmp = load(fullfile(outdir,outstem,seed_name,'results',strcat('final_result_',outstem,'.mat')));
            all_results = tmp.optimal_stats.(metric)';
            seed_alloc = ((seed_num-1)*N_folds + 1) : (seed_num*N_folds);
            acc_vec(:,seed_alloc) = all_results(behav_ind,:);
        end
    case 'LRR' % LRR results
        if exist(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat'))) && store_new == 0
            load(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat')));
        else
            % iterate over behaviours
            for comp = 1:N_behav
                j = behav_ind(comp);
                variable_prefix = strcat('results_var_',num2str(j),'_');
                results_folder = dir(fullfile(outdir,outstem));
                folders = {results_folder.name};
                result_loc = ~cellfun('isempty', strfind(folders,variable_prefix));
                score_dir = folders(result_loc);
                for seed_num = 1:N_seeds
                    seed_name = strcat('seed_', num2str(seed_num));
                    load(fullfile(outdir,outstem,score_dir{1}, seed_name, ...
                        'results','optimal_acc',strcat(outstem,'.mat')));
                    for k = 1:N_folds
                        curr_result_pos = (seed_num-1)*N_folds + k;
                        acc_vec(comp,curr_result_pos) = optimal_statistics{k}.(metric);
                    end
                end
            end
            save(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat')), 'acc_vec')
        end
    case 'Elasticnet' % Elasticnet results
        if exist(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat'))) && store_new == 0
            load(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat')));
        else
            % iterate over behaviours
            for comp = 1:N_behav
                j = behav_ind(comp);
                variable_prefix = strcat('results_var_',num2str(j));
                results_folder = dir(fullfile(outdir,outstem));
                folders = {results_folder.name};
                result_loc = ~cellfun('isempty', strfind(folders,variable_prefix));
                score_dir = folders(result_loc);
                for seed_num = 1:N_seeds
                    seed_name = strcat('seed_', num2str(seed_num));
                    load(fullfile(outdir,outstem,score_dir{1}, seed_name,'optimal_acc', ...
                        strcat(outstem,'_final_acc.mat')));
                    for k = 1:N_folds
                        curr_result_pos = (seed_num-1)*N_folds + k;
                        if ~isnan(optimal_statistics{k}.(metric))
                            acc_vec(comp,curr_result_pos) = optimal_statistics{k}.(metric);
                        else
                            acc_vec(comp,curr_result_pos) = 0;
                        end
                    end
                end
            end
            save(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat')), 'acc_vec');
        end
        
    %% second level models
    case 'multiKRR'
        for seed_num = 1:N_seeds
            seed_name = strcat('seed_', num2str(seed_num));
            tmp = load(fullfile(outdir,outstem,seed_name,'results',strcat('final_result_',outstem,'.mat')));
            multikrr_metric = find(strcmp({tmp.optimal_stats{1,1}.description},metric));
            seed_alloc = ((seed_num-1)*N_folds + 1) : (seed_num*N_folds);
            acc_vec(:,seed_alloc) = tmp.optimal_stats{1,1}(multikrr_metric).value';
        end
        acc_vec = acc_vec(behav_ind,:);
    case 'stacking'
        if exist(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat'))) && store_new == 0
            load(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat')));
        else
            for var = 1:N_behav
                j = behav_ind(var);
                variable = strcat('variable_', num2str(j));
                for seed_num = 1:N_seeds
                    seed_name = strcat('seed_', num2str(seed_num));
                    for outerfold = 1:N_folds
                        outerfold_str = strcat('fold_', num2str(outerfold));
                        tmp = load(fullfile(outdir,outstem,variable, seed_name, outerfold_str, ...
                            'results','optimal_acc', strcat(outstem,'.mat')));
                        acc_vec(var,(seed_num-1)*N_folds + outerfold) = tmp.optimal_statistics{1,1}.(metric);
                    end
                end
            end
            save(fullfile(outdir,outstem,strcat(outstem, '_', metric, '_acc_vec.mat')), 'acc_vec')
        end
        
    %% extract mean and best performing first-level KRR models
    case 'mean' % mean performance over all models of specific modality
        acc_vec = zeros(N_behav,N_seeds*N_folds);
        modality_outstem = select_modality(split_outstem{2});
        for modality = 1:length(modality_outstem)
            modality_name = strcat(split_outstem{3},'_',modality_outstem{modality});
            acc_vec = acc_vec + CBIG_MMP_HCP_read_model_results(modality_name, outdir, ...
                N_seeds, N_folds, behav_ind, metric, store_new)/length(modality_outstem);
        end
    case 'best' % best of 1st level KRR results for each modality
        acc_vec = -ones(N_behav,N_seeds*N_folds);
        modality_outstem = select_modality(split_outstem{2});
        for modality = 1:length(modality_outstem)
            modality_name = strcat(split_outstem{3},'_',modality_outstem{modality});
            new_acc_vec = CBIG_MMP_HCP_read_model_results(modality_name, outdir, ...
                N_seeds, N_folds, behav_ind, metric, store_new);
            % replace if average over first 58 behaviours is higher
            if mean(new_acc_vec(1:58,:), 'all') > mean(acc_vec(1:58,:), 'all')
                acc_vec(1:58,:) = new_acc_vec(1:58,:);
            end
            % replace if factor score is higher
            for factor=59:61
                if mean(new_acc_vec(factor,:)) > mean(acc_vec(factor,:))
                    acc_vec(factor,:) = new_acc_vec(factor,:);
                end
            end
        end
end
end

function modality_outstem = select_modality(selection)
% define all possible features from each modality
% choose features based on selection
switch selection
    case 'struct'
        modality_outstem = {'features_cv' 'features_ca' 'features_ct'};
    case 'tbss'
        modality_outstem = {'features_tbss_FA' 'features_tbss_MD' ...
            'features_tbss_AD' 'features_tbss_RD' 'features_tbss_OD' ...
            'features_tbss_ICVF' 'features_tbss_ISOVF'};
    case 'sc'
        modality_outstem = {'features_schaefer_FA' 'features_schaefer_MD' ...
            'features_schaefer_AD' 'features_schaefer_RD' 'features_schaefer_OD' ...
            'features_schaefer_ICVF' 'features_schaefer_ISOVF' ...
            'features_schaefer_streamcount_log' 'features_schaefer_streamlen' };
    case 'fmri'
        modality_outstem = {'features_rs' 'features_social' 'features_gamb' ...
            'features_lang' 'features_wm' 'features_motor'};
end
end