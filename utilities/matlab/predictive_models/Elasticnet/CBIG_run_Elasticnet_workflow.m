function CBIG_run_Elasticnet_workflow(params)

% This function uses glmnet to perform an elastic-net regression and is
% adapted from CBIG_LRR_workflow_1measure from the CBIG_repo. There are
% 2 parameters that are optimised for elasticnet: lambda (the
% regularisation parameter) and alpha (the ratio of L1 to L2 regularisation). 
% The two parameters are optimized using grid search.
% 
% Inputs:
%   params is a struct containing the following fields:
%   - params.sub_fold (compulsory)
%     A structure containing the information of how the data are separated
%     into training and test sets. 
%     params.sub_fold(i).fold_index is a #subjects x 1 binary vector. 1
%     refers to the corresponding subject is a test subject in the i-th
%     test fold. 0 refers to the corresponding subject is a training
%     subject in the i-th test fold.
% 
%   - params.feature_mat (compulsory)
%     A matrix of the features used as the independent variable in the
%     prediction. It can be a 2-D matrix with dimension of #features x
%     #subjects or a 3-D matrix with dimension of #ROIs1 x #ROIs2 x
%     #subjects. If params.feature_mat is 3-D, it is a connectivity matrix
%     between two sets of ROIs, and only the lower-triangular off-diagonal
%     entries will be considered as features because the connectivity
%     matrix is symmetric.
% 
%   - params.covariates (compulsory)
%     A matrix of the covariates to be regressed out from y (dimension:
%     #subjects x #regressors).
%     If the users only want to demean, an empty matrix should be passed in
%     ; if the users do not want to regress anything, set params.covariates
%     = 'NONE'.
% 
%   - params.y (compulsory)
%     A column vector of the target variable to be predicted (dimension:
%     #subjects x 1).
%
%   - params.outdir (compulsory)
%     A string of the full-path output directory.
%
%   - params.split_name (compulsory)
%     A string of the name of the split in the training test fold.
% 
%   - params.outstem (compulsory)
%     A string appended to the filename to specify the output files (y
%     after regression, accuracy files, ...). For example, if outstem =
%     'AngAffect_Unadj', then the final optimal accuracy file will be names 
%     as
%     fullfile(outdir, 'results', 'optimal_acc', 'AngAffect_Unadj.mat'),
%     and the optimal model hyperparameters of k-th fold will be saved in
%     fullfile(outdir, 'params', ['dis_' k '_cv'], ...
%     'selected_parameters_AngAffect_Unadj.mat').
%
%   - params.glmnet_dir (optional)
%     A string, the full path of the directory to store the glmnet code.
%     Default is
%     $CBIG_CODE_DIR/external_packages/matlab/non_default_packages/glmnet/glmnet_matlab
% 
%   - params.num_inner_folds (optional)
%     A string or a scalar. The number of inner-loop folds within each
%     training fold for hyperparameter selection.
%     If params does not have a field called num_inner_folds, this script
%     uses the same number of outer-loop folds as the number of inner-loop
%     folds.
% 
%   - params.lambda (optional)
%     A vector of regularization weights to search through. The best
%     lambda for each training-test split will be selected from the vector.
%     The default set of vectors is
%     [0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20].
%     Alternatively, glmnet will set lambda for you by passing '[]'.
%
%   - params.alpha (optional)
%     A vector of ratios of L1 to L2 regularization to search through. The
%     best alpha for each training-test split will be selected from the vector.
%     The default set of vectors is
%     [0.01 0.1 0.4 0.7 0.9 0.99]
%
%   - params.metric (optional)
%     A string, specifying how the value of lambda will be optimised with.
%     Can choose the following: corr, COD, predictive_COD, MAE,
%     MAE_norm,MSE,MSE_norm. Metric is predictive_COD by default.
%
%   - params.norm (optional)
%     A logical, specifying if features should be normalized. If true,
%     features will be normalized over subjects in each test/training
%     fold. False by default.
%
%  Outputs:
%   Three files will be generated in the output directory: 'optimal_acc',
%   'params', and 'y'.
%
%   In 'optimal_acc', a mat file with the following fields are generated:   
%   - acc_corr_train
%     Training accuracy (given in correlation).
%
%   - acc_corr_test
%     Test accuracies (given in correlation).
%
%   - y_predict 
%     Predicted target values.
%
%   - y_pred_train
%     Predicted target values of training subjects.
%
%   - optimal_statistics 
%     A cell array of size equal to the number of folds. Each 
%     element in the cell array is a structure storing the accuracies 
%     of each possible accuracy metric (eg. corr, MAE, etc).
%
%   In 'params', a mat file with the following fields are generated: 
%   - selcted_parameters
%     The selected alpha and lambda values for each fold.
%
%   - acc_train
%     The y predictions from each fold of the inner-loop cross validation
%     and the mean accuracy over all inner-loop cross validations.
%
%   In 'y', a mat file with the following fields are generated: 
%   - y_orig
%     Orignal y values passed into the workflow.
%
%   - y_resid
%     Residual y values after regressing out covariates.
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Preparation
if(~isfield(params, 'glmnet_dir'))
    params.glmnet_dir = fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
        'non_default_packages', 'glmnet', 'glmnet_matlab');
end
addpath(params.glmnet_dir)

if(~isfield(params, 'num_inner_folds'))
    params.num_inner_folds = length(params.sub_fold);
end

if(~isfield(params, 'lambda'))
    params.lambda = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
end
lambda_sorted = sort(params.lambda, 'descend'); % sort lambda in descending order for glmnet

if(~isfield(params, 'alpha'))
    params.alpha = [0.01 0.1 0.4 0.7 0.9 0.99];
end

if(~isfield(params, 'metric'))
    params.metric = 'predictive_COD';
end

if(~isfield(params, 'norm'))
    params.norm = false;
end

%% step 1. Regress covariates from y
fprintf('# step 1: regress covariates from y for each fold.\n')
y_dir = fullfile(params.outdir, params.split_name);
CBIG_crossvalid_regress_covariates_from_y( ...
    params.y, params.covariates, params.sub_fold, y_dir, params.outstem);

%% step 2. For each fold, use optimize hyperparameters for elasticnet
fprintf('# step 2: optimize elasticnet regression model and predict target measure.\n')

y_predict = cell(length(params.sub_fold), 1);
y_pred_train = cell(length(params.sub_fold), 1);
acc_corr_test = zeros(length(params.sub_fold), 1);
param_dir = fullfile(params.outdir, params.split_name, '/params');
alpha_set = params.alpha;

for currfold = 1:length(params.sub_fold)
    curr_param_dir = fullfile(param_dir, ['/fold_' num2str(currfold)]);
    if(~exist(curr_param_dir, 'dir'))
        mkdir(curr_param_dir)
    end
    fprintf('=========================== Fold %d ===========================\n', currfold)

    load(fullfile(y_dir, 'y', ['fold_' num2str(currfold)], ['y_regress_' params.outstem '.mat']))
    y_train = y_resid(params.sub_fold(currfold).fold_index==0);
    y_test = y_resid(params.sub_fold(currfold).fold_index==1);

    if(ndims(params.feature_mat) == 3)
        full_feat_train = params.feature_mat(:,:, params.sub_fold(currfold).fold_index==0);
        full_feat_test = params.feature_mat(:,:, params.sub_fold(currfold).fold_index==1);
        full_feat_train = reshape3d_to_2d(full_feat_train);
        full_feat_test = reshape3d_to_2d(full_feat_test);
    else
        full_feat_train = params.feature_mat(:, params.sub_fold(currfold).fold_index==0);
        full_feat_test = params.feature_mat(:, params.sub_fold(currfold).fold_index==1);
    end

    % normalize data if needed
    if params.norm
        full_feat_train = (full_feat_train - mean(full_feat_train,2)) ./ std(full_feat_train,0,2);
        full_feat_test = (full_feat_test - mean(full_feat_test,2)) ./ std(full_feat_test,0,2);
    end 
    
    % no feature selection, so set used features to full features
    feat_train = full_feat_train;
    feat_test = full_feat_test;

    param_file = fullfile(curr_param_dir, ['selected_parameters_' params.outstem '.mat']);
    if(~exist(param_file, 'file'))
    % step 2. select hyperparameters: gridsearch over alpha and lambda over full set of features
        fprintf('>>> Performing gridsearch over alpha and lambda \n') 
        if isempty(lambda_sorted)
            fprintf('>>> Letting Glmnet choose lambda \n') 
        end 
        for alpha_idx = 1:length(alpha_set)
            alpha_train = alpha_set(alpha_idx);
            [acc_metric_tmp{alpha_idx}, neg_loss_tmp{alpha_idx}, new_lambda{alpha_idx},...
        y_pred_tmp{alpha_idx}] = CBIG_Elasticnet_innerloop_cv_glmnet(...
        feat_train', y_train, alpha_train, lambda_sorted, params.num_inner_folds, params.metric);
    end
    best_alpha_idx = find(max(cellfun(@(x)max(x),neg_loss_tmp)),1);
    curr_alpha = alpha_set(best_alpha_idx);
    best_lambda_idx = find(neg_loss_tmp{best_alpha_idx} == max(neg_loss_tmp{best_alpha_idx}),1);
    curr_lambda = new_lambda{best_alpha_idx}(best_lambda_idx);
    save(param_file, 'curr_alpha', 'curr_lambda')
    curr_acc_train = acc_metric_tmp{best_alpha_idx}(best_lambda_idx); 
    curr_y_pred_alpha = y_pred_tmp{best_alpha_idx};
    curr_y_pred = cellfun(@(x)x(:,best_lambda_idx),curr_y_pred_alpha, 'UniformOutput', false);
    save([curr_param_dir '/acc_train_' params.outstem '.mat'], 'curr_acc_train', 'curr_y_pred');
    else
        load(param_file)
        if(~exist([curr_param_dir '/acc_train_' params.outstem '.mat'], 'file'))
            [curr_acc_train,~,~,curr_y_pred] = CBIG_Elasticnet_innerloop_cv_glmnet(...
        feat_train', y_train, curr_alpha, curr_lambda, params.num_inner_folds, params.metric);
        save([curr_param_dir '/acc_train_' params.outstem '.mat'], 'curr_acc_train', 'curr_y_pred');
        else
            load([curr_param_dir '/acc_train_' params.outstem '.mat'])
        end
    end
 
    opt_alpha(currfold, 1) = curr_alpha;
    opt_lambda(currfold, 1) = curr_lambda;
    acc_metric_train(currfold, 1) = curr_acc_train;
    fprintf('>>> Current optimal alpha: %f, lambda: %f, training accuracy: %f\n', ...
        opt_alpha(currfold), opt_lambda(currfold), acc_metric_train(currfold));

    %% step 3. Find outer-loop CV test accuracy
    fprintf('Test fold %d:\n', currfold);

    [acc_corr, y_predict{currfold}, y_pred_train{currfold}, optimal_statistics{currfold}] = ...
        CBIG_Elasticnet_train_test_glmnet(feat_train', feat_test', y_train, y_test, curr_alpha, curr_lambda);
    acc_corr_test(currfold, 1) = acc_corr;

    clear acc_corr
    
end

%% save optimal accuracies per variable
opt_out = fullfile(params.outdir, params.split_name, ...
    'optimal_acc', [params.outstem '_final_acc.mat']);
mkdir(fullfile(params.outdir, params.split_name, 'optimal_acc'))
save(opt_out, 'acc_metric_train', 'acc_corr_test', 'optimal_statistics','y_predict','y_pred_train')
fprintf('Finished!\n')

rmpath(params.glmnet_dir)
end

function out = reshape3d_to_2d(features)
    % reshapes a #ROI x #ROI x #subjects matrix into
    % #ROI x #subjects by extracting the lower triangle
    % of the correlation
    temp = ones(size(features,1), size(features,2));
    tril_ind = tril(temp, -1);
    features_reshaped = reshape(features, size(features,1)*size(features,2), size(features, 3));
    out = features_reshaped(tril_ind==1, :);
end
