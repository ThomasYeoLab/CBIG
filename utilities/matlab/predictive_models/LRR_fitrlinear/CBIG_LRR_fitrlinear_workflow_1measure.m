function [ acc_corr_train, acc_corr_test, y_predict] = CBIG_LRR_fitrlinear_workflow_1measure( params )

% [ acc_corr_train, acc_corr_test, y_predict] = CBIG_LRR_fitrlinear_workflow_1measure( params )
%
% This function is an adaptation to the CBIG_LRR_workflow_1measure in
% CBIG_repo. Instead of the usual inverse of matrices to fit the LRR, this
% function uses the fitrlinear function for fitting which is
% considerably faster. So instead of calling CBIG_LRR_train_test and
% CBIG_LRR_innerloop_cv, an adapted CBIG_LRR_train_test_fitrlinear and
% CBIG_LRR_innerloop_cv_fitrlinear is used.
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
%   - params.num_inner_folds (optional)
%     A string or a scalar. The number of inner-loop folds within each
%     training fold for hyperparameter selection.
%     If params does not have a field called num_inner_folds, this script
%     uses the same number of outer-loop folds as the number of inner-loop
%     folds.
% 
%   - params.outdir (compulsory)
%     A string of the full-path output directory.
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
%
%
%   - params.cov_X (optional)
%     A matrix specifiying the covariates which need to be regressed from features.
%     Each row corresponds to one subject. Each column corresponds to one covariate.
%     If no covariate needs to be regressed from features, either skip this field, 
%     or set params.cov_X = [].
%
%   - params.gpso_dir (optional)
%     A string, the full path of the directory to store the cloned git
%     repository for Gaussian process code. Current script needs the
%     Gaussian-Process Surrogate Optimisation (GPSO;
%     https://github.com/jhadida/gpso) package to optimize the objective
%     function. This git repository and its dependency "deck"
%     (https://github.com/jhadida/deck) need to be both cloned into the
%     same folder, which is passed into current script through
%     params.gpso_dir.
%     Default is
%     $CBIG_CODE_DIR/external_packages/matlab/non_default_packages/Gaussian_Process
% 
%   - params.domain (optional)
%     The searching domain of parameters used by the Gaussian process
%     optimization algorithm (dimension: #hyperparameters x 2). Default is
%     [0.001 0.1; 3 8], where 0.001-0.1 is the searching domain for the
%     feature selection threshold, and 3-8 is the searching domain for the
%     L2 regularization hyperparameter (after taking logarithm). If feature
%     selection is not needed a 1 x 2 vector with the search domain for
%     lambda can be used. For more information, please refer to:
%     https://github.com/jhadida/gpso
% 
%   - params.eval (optional)
%     The maximal evaluation times of the objective function (a scalar).
%     Default is 15. If it is too large, the runtime would be too long; if
%     it is too small, the objective function may not be able to reach its
%     optimal value. The users need to consider this trade-off when setting
%     this variable. For more information, please refer to: 
%     https://github.com/jhadida/gpso
% 
%   - params.tree (optional)
%     A scalar, the depth of the partition tree used to explore
%     hyperparameters. Default is 3. For more information, please refer to:
%     https://github.com/jhadida/gpso
%
%   - params.metric (optional)
%     A string, specifying how the value of lambda will be optimised with.
%     Can choose the following: corr, COD, predictive_COD, MAE,
%     MAE_norm,MSE,MSE_norm. Default is predictive_COD
%
%  output (saved)
%   - acc_corr_train - training accuracy (given in correlation)
%   - acc_corr_test - test accuracies (given in correlation)
%   - y_predict - predicted target values
%   - optimal_statistics - a cell array of size equal to the number of
%   folds. Each element in the cell array is a structure storing the
%   accuracies of each possible accuracy metric (eg. corr, MAE, etc).
% 
% Written by Yanrui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% check if submodules are correctly setup
if(~isfield(params, 'gpso_dir'))
    params.gpso_dir = fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
        'non_default_packages', 'Gaussian_Process');
end
gpso_info = dir(fullfile(params.gpso_dir, 'deck'));
if(isempty(setdiff({gpso_info.name}, {'.', '..'})))
    curr_dir = pwd;
    cd(getenv('CBIG_CODE_DIR'))
    command = sprintf('git submodule update --init --recursive');
    system(command);
    cd(curr_dir)
end

addpath(fullfile(params.gpso_dir, 'gpso'))
addpath(fullfile(params.gpso_dir, 'deck'))
dk_startup();

%% Preparation
if(~isfield(params, 'num_inner_folds'))
    params.num_inner_folds = length(params.sub_fold);
end

if(~isfield(params, 'domain'))
    params.domain = [0.001 0.1; 3 8];
end

if(~isfield(params, 'eval'))
    params.eval = 15;
end

if(~isfield(params, 'tree'))
    params.tree = 3;
end

if(~isfield(params, 'metric'))
    params.metric = 'predictive_COD';
end

%% step 1. Regress covariates from y
fprintf('# step 1: regress covariates from y for each fold.\n')
CBIG_crossvalid_regress_covariates_from_y( ...
    params.y, params.covariates, params.sub_fold, params.outdir, params.outstem);

%% step2. For each fold, use Gaussian process to optimize the linear ridge regression model
param_dir = fullfile(params.outdir, 'params');
fprintf('# step 2: optimize linear ridge regression model and predict target measure.\n')
fprintf('# Evaluations of GP: %d; depth: %d\n', params.eval, params.tree);
y_predict = cell(length(params.sub_fold), 1);
acc_corr_test = zeros(length(params.sub_fold), 1);
for currfold = 1:length(params.sub_fold)
    curr_param_dir = fullfile(param_dir, ['/fold_' num2str(currfold)]);
    if(~exist(curr_param_dir, 'dir'))
        mkdir(curr_param_dir)
    end
    fprintf('=========================== Fold %d ===========================\n', currfold)
    
    load(fullfile(params.outdir, 'y', ['fold_' num2str(currfold)], ['y_regress_' params.outstem '.mat']))
    y_train = y_resid(params.sub_fold(currfold).fold_index==0);
    y_test = y_resid(params.sub_fold(currfold).fold_index==1);
    
    if(ndims(params.feature_mat) == 3)
        feature_train = reshape_3D_features(params.feature_mat(:,:, params.sub_fold(currfold).fold_index==0));
        feature_test = reshape_3D_features(params.feature_mat(:,:, params.sub_fold(currfold).fold_index==1));
    else
        feature_train = params.feature_mat(:, params.sub_fold(currfold).fold_index==0);
        feature_test = params.feature_mat(:, params.sub_fold(currfold).fold_index==1);
    end

    if(isfield(params, 'cov_X') && ~isempty(params.cov_X) && ~strcmpi(params.cov_X, 'none'))
        [feature_train, beta] = CBIG_regress_X_from_y_train(feature_train', ...
            params.cov_X(params.sub_fold(currfold).fold_index==0, :));
        cov_X_train_mean = mean(params.cov_X(params.sub_fold(currfold).fold_index==0, :));
        feature_test = CBIG_regress_X_from_y_test(feature_test', ...
            params.cov_X(params.sub_fold(currfold).fold_index==1, :), beta, cov_X_train_mean);

        feature_train = feature_train';
        feature_test = feature_test';

        save(fullfile(curr_param_dir, 'feature_regress_beta.mat'), 'beta')
    end
    
    %% step 2.1. select hyperparameters: feature selection & inner-loop CV using GPSO
    param_file = fullfile(curr_param_dir, ['selected_parameters_' params.outstem '.mat']);
    if(~exist(param_file, 'file'))
        objfun = @(parameters) select_hyperparameter( feature_train, feature_test, y_train, ...
            parameters, params.num_inner_folds, params.metric );
        obj = GPSO();
        output =  obj.run(objfun, params.domain, params.eval, 'tree', params.tree);
        if size(params.domain,1) == 1
            fprintf('No feature selection\n');
            curr_threshold = 1;
            curr_lambda = 10^output.sol.x(1);
        else
            curr_threshold = output.sol.x(1);
            curr_lambda = 10^output.sol.x(2);
        end
        save(param_file, 'curr_threshold', 'curr_lambda')

        % save accuracy and innerloop y predictions
        if curr_threshold < 1
            [feat_train, ~] = CBIG_FC_FeatSel( feature_train, feature_test, y_train, curr_threshold );
        else
            feat_train = feature_train;
        end
        [curr_acc_train, curr_y_pred] = CBIG_LRR_fitrlinear_innerloop_cv( feat_train', y_train, ...
            curr_lambda, params.num_inner_folds, params.metric );
        save([curr_param_dir '/acc_train_' params.outstem '.mat'], 'curr_acc_train', 'curr_y_pred');
    else
        load(param_file)
        if(~exist([curr_param_dir '/acc_train_' params.outstem '.mat'], 'file'))
            if curr_threshold < 1
                [feat_train, ~] = CBIG_FC_FeatSel( feature_train, feature_test, y_train, curr_threshold );
            else
                if ndims(feature_train) == 3
                    feat_train = reshape_3D_features(feature_train);
                else
                    feat_train = feature_train;
                end
            end
            [curr_acc_train, curr_y_pred] = CBIG_LRR_fitrlinear_innerloop_cv( feat_train', y_train, ...
                curr_lambda, params.num_inner_folds, params.metric );
            save([curr_param_dir '/acc_train_' params.outstem '.mat'], 'curr_acc_train', 'curr_y_pred');
        else
            load([curr_param_dir '/acc_train_' params.outstem '.mat'])
        end
    end
    opt_threshold(currfold, 1) = curr_threshold;
    opt_lambda(currfold, 1) = curr_lambda;
    acc_corr_train(currfold, 1) = curr_acc_train;
    fprintf('>>> Current optimal feature selection threshold: %f, optimal lambda: %f, training accuracy: %f\n', ...
        opt_threshold(currfold), opt_lambda(currfold), acc_corr_train(currfold));
    
    %% Training & test
    fprintf('Test fold %d:\n', currfold);
    
    if curr_threshold == 1
        feat_train = feature_train;
        if ndims(feat_train) == 3
            feat_train = reshape_3D_features(feat_train);
        end
        feat_test = feature_test;
        if ndims(feat_test) == 3
            feat_test = reshape_3D_features(feat_test);
        end
    else
        [feat_train, feat_test] = CBIG_FC_FeatSel( feature_train, feature_test, y_train, curr_threshold );
    end
    [acc_corr, optimal_statistics{currfold},y_predict{currfold}] = CBIG_LRR_fitrlinear_train_test( ...
        feat_train', feat_test', y_train, y_test, curr_lambda, params.metric );
    
    acc_corr_test(currfold, 1) = acc_corr;
    clear acc_corr
end

opt_out = fullfile(params.outdir, 'results', 'optimal_acc', [params.outstem '.mat']);
mkdir(fullfile(params.outdir, 'results', 'optimal_acc'))
save(opt_out, 'acc_corr_train', 'acc_corr_test', 'y_predict', 'optimal_statistics')
fprintf('Finished!\n')

rmpath(fullfile(params.gpso_dir, 'gpso'))
rmpath(fullfile(params.gpso_dir, 'deck'))

end


function out = select_hyperparameter( feature_train, feature_test, y_train, parameters, tot_folds, metric )

if length(parameters) == 1
    FS_threshold = 1;
    lambda = 10^parameters(1);
    feat_train = feature_train;
    if ndims(feat_train) == 3
        feat_train = reshape_3D_features(feat_train);
    end
else
    FS_threshold = parameters(1);
    lambda = 10^parameters(2);
    [feat_train, ~] = CBIG_FC_FeatSel( feature_train, feature_test, y_train, FS_threshold );
end

fprintf('Feature seletion threshold: %f, lambda: %f ...\n', FS_threshold, lambda);

[mean_loss, ~] = CBIG_LRR_fitrlinear_innerloop_cv( feat_train', y_train, lambda, tot_folds, metric );
out = mean_loss;

end

function out = reshape_3D_features(features)
    temp = ones(size(features,1), size(features,2));
    tril_ind = tril(temp, -1);
    
    features_reshaped = reshape(features, size(features,1)*size(features,2), size(features, 3));
    out = features_reshaped(tril_ind==1, :);
end

