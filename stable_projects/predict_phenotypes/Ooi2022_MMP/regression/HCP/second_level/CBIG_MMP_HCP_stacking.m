function CBIG_MMP_HCP_stacking(outerFolds, innerFolds, krr_output_dir, feature_mats_cell, ...
    outstem_name, outdir, split_idx, ymat, covmat, y_idx, second_lvl)

%function CBIG_MMP_HCP_stacking(outerFolds, innerFolds, krr_output_dir, feature_mats_cell, ...
%   outstem_name, outdir, split_idx, ymat, covmat, y_idx, second_lvl)
%
% This function prepares the input parameters for the stacking model using the
% cross validation workflow. User can specify the type of regression to be done for the
% second level.
%
% Inputs:
%   - outerFolds
%     Number of outer-loop cross validation folds.
%
%   - innerFolds
%     Number of inner-loop cross validation folds
%
%   - krr_output_dir
%     Path to a folder containing results of single feature KRR models.
%
%   - feature_mats_cell
%     Cell of base file name for the input features.
%
%   - outstem_name
%     Name of output folder.
%
%   - split_idx
%     The split which the stacking workflow is being run for.
%
%   - outdir
%     The full path of output directory where `outstem_name` will be created.
%
%   - ymat
%     Name of the output mat file containing behavioural variables for all
%     subjects in the regression analysis.
%
%   - covmat
%     Name of the output mat file containing covariates for all subjects in the
%     regression analysis.
%
%   - y_idx
%     The behaviour index which the stacking workflow is being run for.
%
%   - second_lvl
%     Type of regression model to run for stacking. Choices are between
%     'KRR', 'LRR' and 'Elasticnet'.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% format params that are read in
% folds
seed = split_idx;
seed_name = strcat('seed_', num2str(seed));
y_name = strcat('variable_', num2str(y_idx));
num_outer_folds = outerFolds;
param.num_inner_folds = innerFolds;
% feature details
outstem = outstem_name;
param.outstem = outstem;
param.split_name = seed_name;
% subject details
y_mat = ymat;
cov_mat = covmat;
% regression settings
lvl1_lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
lambda_set = [ 0 ];
lambda_domain = [-8 -4];
param.lambda = lambda_set;
param.lambda_set = lambda_set;
param.domain = lambda_domain;
param.with_bias = 1;
param.ker_param.type = 'corr';
param.ker_param.scale = NaN;
param.threshold_set = [];
param.metric = 'predictive_COD';

for fold_idx = 1:num_outer_folds
    % set up fold specific variables
    fold_name = strcat('fold_', num2str(fold_idx));
    param.outdir = fullfile(outdir, outstem, y_name, seed_name, fold_name);
    fprintf('--- Outerfold %d... ---\n', fold_idx)
    %% get subfold
    fprintf('[1] Fetch subfold... \n')
    feature_name = feature_mats_cell{1};
    feature_outstem = strcat('KRR_',feature_name);
    subfold = load(fullfile(krr_output_dir,feature_outstem, seed_name,'results','no_relative_10_fold_sub_list.mat'));
    param.sub_fold = subfold.sub_fold(fold_idx);
    
    %% generate y matrix
    fprintf('[2] Fetch y matrix... \n')
    fprintf('Using existing y file \n')
    y_temp = load(fullfile(krr_output_dir,y_mat));
    y = y_temp.y;
    param.y = y(:,y_idx);
    
    %% generate covariate matrix
    fprintf('[3] Fetch covariate matrix... \n')
    fprintf('Using existing covariate file \n')
    cov_temp = load(fullfile(krr_output_dir, cov_mat));
    cov = cov_temp.covariates;
    param.covariates = cov;
    
    %% generate features
    fprintf('[4] Loading 1st level predictions.. \n')
    
    all_feat_preds = [];
    for n = 1:size(feature_mats_cell,2)
        feature_name = feature_mats_cell{n};
        fprintf('%s \n', feature_name);
        feature_outstem = strcat('KRR_',feature_name);
        feature_result_mat = strcat('acc_KRR_',feature_name,'.mat');
        feat_train_tmp = load(fullfile(krr_output_dir,feature_outstem,seed_name, ...
            'results','innerloop_cv',fold_name,feature_result_mat));
        feat_test = load(fullfile(krr_output_dir,feature_outstem,seed_name, ...
            'results','test_cv',fold_name,feature_result_mat));
        feat_opt = load(fullfile(krr_output_dir,feature_outstem,seed_name, ...
            'results', strcat('final_result_',feature_outstem,'.mat')));
        feat_lambda = feat_opt.optimal_lambda(fold_idx, y_idx);
        % unscramble feat_train - needs to be updated if KRR cv rng changes
        rng(1, 'twister');
        cv_idx = cvpartition(sum(~param.sub_fold.fold_index), 'kfold', param.num_inner_folds);
        scrambled_idx = 1;
        clear feat_train
        for i = 1:param.num_inner_folds
            last_idx = scrambled_idx + cv_idx.TestSize(i) - 1;
            feat_train(cv_idx.test(i)) = ...
                feat_train_tmp.y_pred{:,find(lvl1_lambda_set == feat_lambda)}{:,y_idx}(scrambled_idx:last_idx);
            scrambled_idx = last_idx + 1;
        end
        feat_pred_y(~param.sub_fold.fold_index) = feat_train;
        feat_pred_y(param.sub_fold.fold_index) = feat_test.y_p{:,find(lvl1_lambda_set == feat_lambda)}{:,y_idx};
        all_feat_preds = [all_feat_preds; feat_pred_y];
    end
    param.feature_mat = all_feat_preds;
    
    %% choose regression workflow
    if strcmp(second_lvl,'KRR')
        fprintf('[5] Run KRR workflow...\n')
        CBIG_KRR_workflow(param);
    elseif strcmp(second_lvl,'LRR')
        fprintf('[5] Run LRR workflow...\n')
        CBIG_LRR_fitrlinear_workflow_1measure(param);
    elseif strcmp(second_lvl,'elasticnet')
        fprintf('[5] Run elasticnet workflow...\n')
        CBIG_run_Elasticnet_workflow(param);
    end
end
