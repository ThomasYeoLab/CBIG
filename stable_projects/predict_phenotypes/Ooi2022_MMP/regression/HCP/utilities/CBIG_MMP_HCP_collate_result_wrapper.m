function CBIG_MMP_HCP_collate_result_wrapper(results_dir, save_dir, reg, metric)

% CBIG_MMP_ABCD_generate_collate_wrapper(results_dir, save_dir, reg, metric)
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
behav_ind = [1:61];
varimax_cog = 60;
varimax_dis = 59;
varimax_emo = 61;
all_scores = 1:58;
N_seed = 60;
N_fold = 10;
groups = {'t1' 'tbss' 'tractography' 'fmri'};
N_groups = length(groups);
store_new = 0;

% create output directory
if ~exist(save_dir)
    mkdir(save_dir)
end

%% base features
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

%% define correct outstem for each group based on type of regression and read results
% read single feature model results
if ~contains(reg,'models')
    % set scores and outstems
    score_name = {'Cognition', 'Dissatisfaction', 'Emotional Rec' 'Avg'};
    results.score_name = score_name;
    N_score = length(score_name);
    best_models = cell(N_groups,N_score);
    score_fold_task_modalitymean = zeros(length(behav_ind),N_seed,N_groups);
    score_fold_task_best = zeros(N_score,N_seed,N_groups);
    for n = 1:N_groups
        switch n
            case 1
                outstem = all_outstems.t1;
            case 2
                outstem = all_outstems.tbss;
            case 3
                outstem = all_outstems.tractography;
            case 4
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
            acc_vec_tmp = CBIG_MMP_HCP_read_model_results(outstem{i}, results_dir, ...
                N_seed, N_fold, behav_ind, metric, store_new);
            acc_score_fold_task(:,:,i) = mean_over_seed(acc_vec_tmp, N_seed, N_fold);
        end
        
        score_fold_task = zeros(N_score,N_seed,N_task);
        score_fold_task(1,:,:) = squeeze(mean(acc_score_fold_task(varimax_cog,:,:),1));
        score_fold_task(2,:,:) = squeeze(mean(acc_score_fold_task(varimax_dis,:,:),1));
        score_fold_task(3,:,:) = squeeze(mean(acc_score_fold_task(varimax_emo,:,:),1));
        score_fold_task(4,:,:) = squeeze(mean(acc_score_fold_task(all_scores,:,:),1));
        
        % save to struct
        results.(groups{n}) = score_fold_task;
        % save mean scores
        score_fold_task_modalitymean(:,:,n) = mean(acc_score_fold_task,3);
        % save best scores
        best_score_loc = squeeze(mean(score_fold_task,2)) == max(squeeze(mean(score_fold_task,2)),[], 2);
        for g = 1:N_score
            best_models{n,g} = outstem{find(best_score_loc(g,:))};
            score_fold_task_best(g,:,n) = score_fold_task(g,:,find(best_score_loc(g,:)));
        end
    end
    
    % save outstems
    results.all_outstems = all_outstems;
    % save best scores to struct
    results.best = score_fold_task_best;
    results.best_models = best_models;
    % save mean scores to struct, and calculate mean over each modality
    results.mean_all = score_fold_task_modalitymean;
    N_task = 4;
    score_fold_task = zeros(N_score,N_seed,N_task);
    score_fold_task(1,:,:) = squeeze(mean(score_fold_task_modalitymean(varimax_cog,:,:),1));
    score_fold_task(2,:,:) = squeeze(mean(score_fold_task_modalitymean(varimax_dis,:,:),1));
    score_fold_task(3,:,:) = squeeze(mean(score_fold_task_modalitymean(varimax_emo,:,:),1));
    score_fold_task(4,:,:) = squeeze(mean(score_fold_task_modalitymean(all_scores,:,:),1));
    results.mean = score_fold_task;
    save(fullfile(save_dir,strcat(reg,'_', metric,'_results.mat')), 'results');
    
% read combined feature model results
elseif strcmp(reg, 'combined_models')
    % set scores and outstems
    score_name = {'Cognition', 'Dissatisfaction', 'Emotional Rec'};
    results.score_name = score_name;
    N_score = length(score_name);
    outstem = {'best' 'multiKRR_k1_rs_k2_wm_k3_lang_k4_gamb_k5_social_k6_motor'...
        'stacking_LRR_rs_wm_lang_gamb_social_motor' 'stacking_LRR_all'};
    N_task = length(outstem);
    
    % read results
    acc_score_fold_task = zeros(length(behav_ind),N_seed,N_task-1);
    for i = 2:N_task
        acc_vec_tmp = CBIG_MMP_HCP_read_model_results(outstem{i}, results_dir, ...
            N_seed, N_fold, behav_ind, metric, store_new);
        acc_score_fold_task(:,:,i-1) = mean_over_seed(acc_vec_tmp, N_seed, N_fold);
    end
    
    score_fold_task = zeros(N_score,N_seed,N_task);
    score_fold_task(1,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_cog,:,:),1));
    score_fold_task(2,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_dis,:,:),1));
    score_fold_task(3,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_emo,:,:),1));
    
    % read best modality from KRR mat file
    indiv_res = load(fullfile(save_dir, strcat('KRR_',metric,'_results.mat')));
    best_score_loc = squeeze(mean(indiv_res.results.best,2)) ...
        == max(squeeze(mean(indiv_res.results.best,2)),[], 2);
    for g = 1:N_score
        score_fold_task(g,:,1) = indiv_res.results.best(g,:,find(best_score_loc(g,:)));
    end
    
    % save outstems
    results.all_outstems = outstem;
    % save results to mat file
    results.combined = score_fold_task;
    save(fullfile(save_dir,strcat(reg,'_',metric,'_results.mat')), 'results');
    
elseif strcmp(reg, 'best_models')
    % set scores and outstems
    score_name = {'Cognition', 'Dissatisfaction', 'Emotional Rec'};
    results.score_name = score_name;
    N_score = length(score_name);
    outstem = {'best' 'stacking_LRR_rs_wm_lang_gamb_social_motor' ...
        'stacking_LRR_best_cog' 'stacking_LRR_best_satisf' ...
        'stacking_LRR_best_er' 'stacking_LRR_all'};
    N_task = length(outstem);
    
    % read results
    acc_score_fold_task = zeros(length(behav_ind),N_seed,N_task-1);
    for i = 2:N_task
        acc_vec_tmp = CBIG_MMP_HCP_read_model_results(outstem{i}, results_dir, ...
            N_seed, N_fold, behav_ind, metric, store_new);
        acc_score_fold_task(:,:,i-1) = mean_over_seed(acc_vec_tmp, N_seed, N_fold);
    end
    
    score_fold_task = zeros(N_score,N_seed,N_task);
    score_fold_task(1,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_cog,:,:),1));
    score_fold_task(2,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_dis,:,:),1));
    score_fold_task(3,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_emo,:,:),1));
    
    % read best modality from KRR mat file
    indiv_res = load(fullfile(save_dir, strcat('KRR_',metric,'_results.mat')));
    best_score_loc = squeeze(mean(indiv_res.results.best,2)) ...
        == max(squeeze(mean(indiv_res.results.best,2)),[], 2);
    for g = 1:N_score
        score_fold_task(g,:,1) = indiv_res.results.best(g,:,find(best_score_loc(g,:)));
    end
    
    % save outstems
    results.all_outstems = outstem;
    % save results to mat file
    results.best = score_fold_task;
    save(fullfile(save_dir,strcat(reg,'_',metric,'_results.mat')), 'results');
end

end

function [acc_vec_mean] = mean_over_seed(acc_vec_tmp, N_seed, N_fold)
% get the mean over outer folds for each seed
for seed = 1:N_seed
    acc_vec_mean(:,seed) = mean(acc_vec_tmp(:,[(seed-1)*N_fold+1:seed*N_fold]),2);
end
end
