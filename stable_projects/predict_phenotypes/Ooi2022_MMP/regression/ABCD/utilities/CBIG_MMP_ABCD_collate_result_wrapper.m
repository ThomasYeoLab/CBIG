function CBIG_MMP_ABCD_collate_result_wrapper(results_dir, save_dir, reg, metric)

% CBIG_MMP_ABCD_generate_collate_wrapper(results_dir, save_dir, reg, metric)
%
% This function reads the results of KRR, LRR, Elasticnet, multiKRR and stacking. The result is
% saved in a matrix of #behav x #splits. Results are saved for stacking since
% computational time for collating results from these models are high.
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
behav_ind = [1:39];
varimax_cog_ind = 37;
varimax_per_ind = 39;
varimax_psy_ind = 38;
all_scores = 1:36;
N_fold = 120;
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
    'features_tbss_AD' 'features_tbss_RD' ...
    'features_tbss_OD' 'features_tbss_ICVF' ...
    'features_tbss_ISOVF'};
all_outstems.tractography = {'features_schaefer_FA' 'features_schaefer_MD' 'features_schaefer_AD' ...
    'features_schaefer_RD' 'features_schaefer_OD' 'features_schaefer_ICVF' ...
    'features_schaefer_ISOVF' 'features_schaefer_streamcount_log' ...
    'features_schaefer_streamlen'};
all_outstems.fmri = {'features_rs' 'features_nback' 'features_mid' 'features_sst'};

%% define correct outstem for each group based on type of regression and read results
% read single feature model results
if ~contains(reg,'models')
    % set scores and outstems
    score_name = {'Cognition', 'Personality', 'Mental Health' 'Avg'};
    results.score_name = score_name;
    N_score = length(score_name);
    best_models = cell(N_groups,N_score);
    score_fold_task_modalitymean = zeros(length(behav_ind),N_fold,N_groups);
    score_fold_task_best = zeros(N_score,N_fold,N_groups);
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
        acc_score_fold_task = zeros(length(behav_ind),N_fold,N_task);
        for i = 1:length(outstem)
            acc_score_fold_task(:,:,i) = CBIG_MMP_ABCD_read_model_results(outstem{i}, ....
                results_dir, N_fold, behav_ind, metric, store_new);
        end
        
        score_fold_task = zeros(length(score_name),N_fold,N_task);
        score_fold_task(1,:,:) = squeeze(mean(acc_score_fold_task(varimax_cog_ind,:,:),1));
        score_fold_task(2,:,:) = squeeze(mean(acc_score_fold_task(varimax_per_ind,:,:),1));
        score_fold_task(3,:,:) = squeeze(mean(acc_score_fold_task(varimax_psy_ind,:,:),1));
        score_fold_task(4,:,:) = squeeze(mean(acc_score_fold_task(all_scores,:,:),1));
        
        % save individual scores to struct
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
    score_fold_task = zeros(length(score_name),N_fold,N_task);
    score_fold_task(1,:,:) = squeeze(mean(score_fold_task_modalitymean(varimax_cog_ind,:,:),1));
    score_fold_task(2,:,:) = squeeze(mean(score_fold_task_modalitymean(varimax_per_ind,:,:),1));
    score_fold_task(3,:,:) = squeeze(mean(score_fold_task_modalitymean(varimax_psy_ind,:,:),1));
    score_fold_task(4,:,:) = squeeze(mean(score_fold_task_modalitymean(all_scores,:,:),1));
    results.mean = score_fold_task;
    save(fullfile(save_dir,strcat(reg,'_', metric,'_results.mat')), 'results');
    
% read combined feature model results
elseif strcmp(reg, 'combined_models')
    % set scores and outstems
    score_name = {'Cognition', 'Personality', 'Mental Health'};
    results.score_name = score_name;
    N_score = length(score_name);
    outstem = {'best' 'multiKRR_k1_rs_k2_mid_k3_sst_k4_nback' ...
        'stacking_LRR_rs_mid_sst_nback' 'stacking_LRR_all'};
    N_task = length(outstem);
    
    % read results
    acc_score_fold_task = zeros(length(behav_ind),N_fold,N_task-1);
    for i = 2:N_task % save best tasks later
        acc_score_fold_task(:,:,i-1) = CBIG_MMP_ABCD_read_model_results(outstem{i}, ...
            results_dir, N_fold, behav_ind, metric, store_new);
    end
    
    score_fold_task = zeros(N_score,N_fold,N_task);
    score_fold_task(1,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_cog_ind,:,:),1));
    score_fold_task(2,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_per_ind,:,:),1));
    score_fold_task(3,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_psy_ind,:,:),1));
    
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
    
% read combined feature model results
elseif strcmp(reg, 'best_models')
    % set scores and outstems
    score_name = {'Cognition', 'Personality', 'Mental Health'};
    results.score_name = score_name;
    N_score = length(score_name);
    outstem = {'best' 'stacking_LRR_rs_mid_sst_nback' ...
        'stacking_LRR_best_cog' 'stacking_LRR_best_pers' ...
        'stacking_LRR_best_mental' 'stacking_LRR_all'};
    N_task = length(outstem);
    
    % read results
    acc_score_fold_task = zeros(length(behav_ind),N_fold,N_task-1);
    for i = 2:N_task % save best tasks later
        acc_score_fold_task(:,:,i-1) = CBIG_MMP_ABCD_read_model_results(outstem{i}, ...
            results_dir, N_fold, behav_ind, metric, store_new);
    end
    
    score_fold_task = zeros(N_score,N_fold,N_task);
    score_fold_task(1,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_cog_ind,:,:),1));
    score_fold_task(2,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_per_ind,:,:),1));
    score_fold_task(3,:,2:N_task) = squeeze(mean(acc_score_fold_task(varimax_psy_ind,:,:),1));
    
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

