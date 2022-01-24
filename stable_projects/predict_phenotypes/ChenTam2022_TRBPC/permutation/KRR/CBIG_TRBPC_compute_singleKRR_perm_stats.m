function CBIG_TRBPC_compute_singleKRR_perm_stats(singleKRR_dir, outstem, score_ind, perm_seed_start, ...
N_perm, site_list, outdir)

% CBIG_TRBPC_compute_singleKRR_perm_stats(singleKRR_dir,outstem,score_ind,perm_seed_start,N_perm,site_list,outdir)
%
% This function performs the permutation test for the singleKRR model and
% saves out the prediction accuracies of the permutation test
%
% Inputs:
%   - singleKRR_dir:
%     The directory wherre the singleKRR results are
%
%   - outstem
%     A string appended to the filename to specify the output files in the
%     singleKRR prediction. Must be the same as the input argument `outstem`
%     when you perform singleKRR
%
%   - score_ind
%     A scalar. The index of the score you want to perfom the permutation test.
%     Range from 1 to # Target Variables in your kernel regression
%
%   - perm_seed_start
%     A scalar. The function will do the permutation N_perm times using
%     permutation seed from perm_seed_start to perm_seed_start+N_perm-1
%
%   - N_perm
%     Number of permutation
%
%   - site_list
%     A text file where each line is the site ID of one subject. Permutation
%     will only be performed within sites
%
%   - outdir
%     Output directory for the permutation results
%
% Outputs:
%    prediction accuracies of the permutation test will be saved in the
%    outdir
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% check if output already exist
result_file = fullfile(outdir,['/acc_score' num2str(score_ind) '_allFolds_permStart' num2str(perm_seed_start) '.mat']);
if ~exist(result_file,'file')
    %% check input arguments
    if(isstr(score_ind))
        score_ind = str2num(score_ind);
    end
    
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    %% load site ID 
    site_all = CBIG_text2cell(site_list);
    [~, ~, site_ind] = unique(site_all);
    N_site = max(site_ind);
    
    %% load cross-validation split
    sub_fold_file = dir(fullfile(singleKRR_dir,'*sub_list.mat'));
    if length(sub_fold_file) ~= 1
        error('There should be one and only one sub_fold file in the singleKRR directory');
    end
    load(fullfile(singleKRR_dir,sub_fold_file.name));
    N_fold = length(sub_fold);
    
    %% load FSM
    load(fullfile(singleKRR_dir,'FSM','FSM_corr.mat'));
    
    %% compute permutation
    metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
    for i = 1:length(metrics)
        stats_perm.(metrics{i}) = zeros(N_fold,N_perm);
    end
        
    load(fullfile(singleKRR_dir, ['covariates_' outstem '.mat']));
    load(fullfile(singleKRR_dir, ['final_result_' outstem '.mat']));
    for i = 1:N_fold
        disp(['fold ' num2str(i)]);
        % load data
        load([singleKRR_dir '/y/fold_' num2str(i) '/y_regress_all_score.mat']);
        y_resid = y_resid(:,score_ind);
        opt_lambda = optimal_lambda(i,score_ind);
        test_ind = sub_fold(i).fold_index;
        train_ind = ~test_ind;
        
        FSM_train = FSM(train_ind,train_ind,:);
        FSM_pred = FSM(test_ind,train_ind,:);
        
        % number of subjects
        num_train_subj = size(FSM_train,1);
        num_test_subj = size(FSM_pred,1);
        
        % Perform training and prediction
        K_r = FSM_train + opt_lambda*eye(num_train_subj);
        % training
        X = ones(num_train_subj,1);
        inv_K_r = inv(K_r);
        beta_pre = (X' * (inv_K_r * X)) \ X' * inv_K_r;
        for j = 1:N_perm
            % permute behavior within site
            rng(j+perm_seed_start-1);
            y_perm = y_resid;
            for k = 1:N_site
                y_tmp = y_perm(site_ind == k);
                y_tmp = y_tmp(randperm(length(y_tmp)));
                y_perm(site_ind == k) = y_tmp;
            end
            y_train_resid = y_perm(train_ind);
            y_test_resid = y_perm(test_ind);
            % predict
            beta = beta_pre * y_train_resid;
            alpha = inv_K_r * (y_train_resid - X * beta);
            
            y_predicted = FSM_pred * alpha + ones(num_test_subj,1) .* beta;
            for k = 1:length(metrics)
                stats_perm.(metrics{k})(i,j) = ...
                    CBIG_compute_prediction_acc_and_loss(y_predicted,y_test_resid,metrics{k},y_train_resid);
            end
        end
    end
    save(fullfile(outdir,...
        ['/acc_score' num2str(score_ind) '_allFolds_permStart' num2str(perm_seed_start) '.mat']),'stats_perm');
end
end
