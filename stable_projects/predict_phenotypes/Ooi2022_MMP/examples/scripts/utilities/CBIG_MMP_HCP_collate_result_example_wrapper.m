function CBIG_MMP_HCP_collate_result_example_wrapper(results_dir, save_dir, reg, metric)

% CBIG_MMP_HCP_collate_result_example_wrapper(results_dir, save_dir, reg, metric)
%
% This function reads the results of KRR, LRR, Elasticnet, multiKRR and stacking. The result is
% saved in a matrix of #behav x #splits. Results are saved for LRR, Elasticnet and stacking since
% computational time for these are high.
%
% Inputs:
%   - results_dir
%     Directory where results to be read are saved.
%
%   - save_dir
%     Directory to save mat files of results.
%
%   - reg
%     Regression model to be read. Can be chosen from
%     {'KRR', 'LRR', 'Elasticnet', 'combined_models' 'best_models'}.
%
%   - metric
%     Metric to be read. Can be chosen from
%     {'corr', 'COD', 'predictive_COD', 'MAE' 'MAE_norm', 'MSE', 'MSE_norm'}.
%
% Outputs:
%   - acc_vec
%     A matrix of accuracies in the dimensions of #behav x #splits.
%
% Written by Leon_Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% basic fields
behav_ind = [1:3];
varimax_cog = 2;
varimax_dis = 1;
varimax_emo = 3;
N_seed = 1;
N_fold = 3;
N_groups = {'t1' 'fmri'};
store_new = 0;

% create output directory
if ~exist(save_dir)
    mkdir(save_dir)
end

%% base features
all_outstems.t1 = {'features_cv'};
all_outstems.fmri = {'features_rs'};

%% define correct outstem for each group based on type of regression and read results
% read single feature model results
if ~contains(reg,'models')
    % set scores and outstems
    score_name = {'Cognition', 'Dissatisfaction', 'Emotional Rec'};
    results.score_name = score_name;
    N_score = length(score_name);
    
    for n = 1:length(N_groups)
        switch n
            case 1
                outstem = all_outstems.t1;
            case 2
                outstem = all_outstems.fmri;
        end
        
        % append reg as prefix to oustems
        for i = 1:length(outstem)
            outstem{i} = strcat(reg,'_',outstem{i});
        end
        N_task = length(outstem);
        
        % read results
        acc_score_fold_task = zeros(length(behav_ind),N_seed,N_task);
        for i = 1:N_task
            acc_vec_tmp = CBIG_MMP_HCP_read_model_results_example(outstem{i}, results_dir, ...
                N_seed, N_fold, behav_ind, metric, store_new);
            acc_score_fold_task(:,:,i) = mean_over_seed(acc_vec_tmp, N_seed, N_fold);
        end
        
        score_fold_task = zeros(N_score,N_seed,N_task);
        score_fold_task(1,:,:) = squeeze(mean(acc_score_fold_task(varimax_cog,:,:),1));
        score_fold_task(2,:,:) = squeeze(mean(acc_score_fold_task(varimax_dis,:,:),1));
        score_fold_task(3,:,:) = squeeze(mean(acc_score_fold_task(varimax_emo,:,:),1));
        
        % save to struct
        results.(N_groups{n}) = score_fold_task;
    end
    
    % save outstems
    results.all_outstems = all_outstems;
    save(fullfile(save_dir,strcat(reg,'_', metric,'_results.mat')), 'results');
    
% read combined feature model results
elseif strcmp(reg, 'combined_models')
    % set scores and outstems
    score_name = {'Cognition', 'Dissatisfaction', 'Emotional Rec'};
    results.score_name = score_name;
    N_score = length(score_name);
    outstem = {'KRR_features_rs' 'multiKRR_k1_cv_k2_rs'...
        'stacking_LRR_cv_rs'};
    N_task = length(outstem);
    
    % read results
    acc_score_fold_task = zeros(length(behav_ind),N_seed,N_task);
    for i = 1:N_task
        acc_vec_tmp = CBIG_MMP_HCP_read_model_results_example(outstem{i}, results_dir, ...
            N_seed, N_fold, behav_ind, metric, store_new);
        acc_score_fold_task(:,:,i) = mean_over_seed(acc_vec_tmp, N_seed, N_fold);
    end
    
    score_fold_task = zeros(N_score,N_seed,N_task);
    score_fold_task(1,:,:) = squeeze(mean(acc_score_fold_task(varimax_cog,:,:),1));
    score_fold_task(2,:,:) = squeeze(mean(acc_score_fold_task(varimax_dis,:,:),1));
    score_fold_task(3,:,:) = squeeze(mean(acc_score_fold_task(varimax_emo,:,:),1));
    
    % save outstems
    results.all_outstems = outstem;
    % save results to mat file
    results.combined = score_fold_task;
    save(fullfile(save_dir,strcat(reg,'_',metric,'_results.mat')), 'results');
end

end

function [acc_vec_mean] = mean_over_seed(acc_vec_tmp, N_seed, N_fold)
% get the mean over outer folds for each seed
for seed = 1:N_seed
    acc_vec_mean(:,seed) = mean(acc_vec_tmp(:,[(seed-1)*N_fold+1:seed*N_fold]),2);
end
end
