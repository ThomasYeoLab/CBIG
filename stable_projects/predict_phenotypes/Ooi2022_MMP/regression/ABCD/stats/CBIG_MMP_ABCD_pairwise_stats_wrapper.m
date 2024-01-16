function [pairwise_stats, p_single_models, p_combined_models] = CBIG_MMP_ABCD_pairwise_stats_wrapper(outdir)

% [pairwise_stats, p_single_models, p_combined_models] = CBIG_MMP_ABCD_pairwise_stats_wrapper
%
% Wrapper function that calculates the p-value of whether a model performs better than another .
% Assumes that all regression outputs are in the same directory.
%
% Inputs:
%
%   - outdir
%     Directory to results of regression models.
%
% Outputs:
%   - pairwise_stats
%     A struct containing the pairwise p-values for each combination of models.
%
%   - p_single_models
%     A matrix of #N_behav x #models. P-value for each model indicating in predicting
%     specified behaviour. Currently set to calculate p-values for each factor score,
%     and the grand average of all original behaviour prediction results, meaning
%     4 p-values per model is produced.
%
%   - p_combined_models
%     A matrix of #N_behav x #models. P-value for each model indicating in predicting
%     specified behaviour. Currently set to calculate p-values for each factor score.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set utility directory
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'))
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'utilities'))

%% Mean and best of modality comparisons
% list of models
models = {'mean_struct_KRR' 'mean_tbss_KRR' 'mean_sc_KRR' 'mean_fmri_KRR' ...
    'mean_struct_LRR' 'mean_tbss_LRR' 'mean_sc_LRR' 'mean_fmri_LRR' ...
    'mean_struct_Elasticnet' 'mean_tbss_Elasticnet' 'mean_sc_Elasticnet' 'mean_fmri_Elasticnet' ...
    'best_struct_KRR' 'best_tbss_KRR' 'best_sc_KRR' 'best_fmri_KRR'};
metric = 'corr';
N_folds = 120;
behav_ind = [1:39];
N_behav = 4;
store_new = 0;

%% Extract model accuracies and calculate p-value
clear p
m1_acc_vecs = zeros(N_behav,N_folds);
m2_acc_vecs = zeros(N_behav,N_folds);
for i = 1:length(models)
    % get first model accuracies
    outstem_1 = models{i};
    m1_acc_vecs_tmp = CBIG_MMP_ABCD_read_model_results(outstem_1, outdir, N_folds, behav_ind, metric, store_new);
    % reshape into 3 factor scores, and average accuracies for original behaviours
    m1_acc_vecs(1:3,:) = m1_acc_vecs_tmp(37:39,:);
    m1_acc_vecs(4,:) = mean(m1_acc_vecs_tmp(1:36,:));
    
    % get second model accuracies
    for j = 1:length(models)
        if j ~= i
            outstem_2 = models{j};
            m2_acc_vecs_tmp = CBIG_MMP_ABCD_read_model_results(outstem_2, outdir, ...
                N_folds, behav_ind, metric, store_new);
            % reshape into 3 factor scores, and average accuracies for original behaviours
            m2_acc_vecs(1:3,:) = m2_acc_vecs_tmp(37:39,:);
            m2_acc_vecs(4,:) = mean(m2_acc_vecs_tmp(1:36,:));
            
            % get p value
            for k = 1:N_behav
                diff_vec = (m1_acc_vecs(k,:) - m2_acc_vecs(k,:))';
                p(j,k) = CBIG_corrected_resampled_ttest(diff_vec, 3/7, 0);
            end
        else
            p(j,:) = NaN(1,N_behav);
        end
    end
    pairwise_stats.(models{i}) = p;
end

%% Combined feature model comparisons
% list of models
models = {'best_fmri_KRR' 'multiKRR_k1_rs_k2_mid_k3_sst_k4_nback'...
    'stacking_LRR_rs_mid_sst_nback' 'stacking_LRR_all'};
all_metrics = {'corr' 'COD'};
N_behav = 3;

%% Extract model accuracies and calculate p-value
clear p
m1_acc_vecs = zeros(N_behav,N_folds);
m2_acc_vecs = zeros(N_behav,N_folds);
for n = 1:length(all_metrics)
    metric = all_metrics{n};
    for i = 1:length(models)
        % get first model accuracies
        outstem_1 = models{i};
        m1_acc_vecs_tmp = CBIG_MMP_ABCD_read_model_results(outstem_1, outdir, N_folds, behav_ind, metric, store_new);
        m1_acc_vecs(1:3,:) = m1_acc_vecs_tmp(37:39,:);
        
        % get second model accuracies
        for j = 1:length(models)
            if j ~= i
                outstem_2 = models{j};
                m2_acc_vecs_tmp = ...
                    CBIG_MMP_ABCD_read_model_results(outstem_2, outdir, N_folds, behav_ind, metric, store_new);
                m2_acc_vecs(1:3,:) = m2_acc_vecs_tmp(37:39,:);
                
                % get p value
                for k = 1:N_behav
                    diff_vec = (m1_acc_vecs(k,:) - m2_acc_vecs(k,:))';
                    p(j,k) = CBIG_corrected_resampled_ttest(diff_vec, 3/7, 0);
                end
            else
                p(j,:) = NaN(1,N_behav);
            end
        end
        pairwise_stats.(strcat(models{i},'_',metric)) = p;
    end
end

% extract pairwise stats that are relevant
p_single_models = [pairwise_stats.mean_fmri_KRR(1:3,:); pairwise_stats.mean_fmri_LRR(5:7,:); ...
    pairwise_stats.mean_fmri_Elasticnet(9:11,:); pairwise_stats.best_fmri_KRR(13:15,:)];
p_combined_models = [pairwise_stats.best_fmri_KRR_corr(2:4,:); ...
    pairwise_stats.stacking_LRR_rs_mid_sst_nback_corr([2 4],:); ...
    pairwise_stats.best_fmri_KRR_corr(2:4,:); pairwise_stats.stacking_LRR_rs_mid_sst_nback_COD([2 4],:)];

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'))
rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'utilities'))
