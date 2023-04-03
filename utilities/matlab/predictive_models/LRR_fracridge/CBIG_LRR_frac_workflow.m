function CBIG_LRR_frac_workflow( params )

% [ acc_corr_train, acc_corr_test, y_predict] = CBIG_LRR_frac_workflow( params )
%
% This function uses fractional ridge regression, which is a faster version
% of normal LRR. This script allows user to perform prediction across multiple
% behavioral measures in parallel. The original LRR code requires the user to
% pass in each single behaivoral measure at a time, which is time comsuming.
% This script follows the same workflow as the original LRR code LinearRidgeRegression
% and LRR_fitrliear. The only difference is that this script uses a different
% hyperparameter.
%
% For more information about fractional LRR, please refer to:
% Ariel Rokem, Kendrick Kay, Fractional ridge regression: a fast, interpretable
% reparameterization of ridge regression, GigaScience, Volume 9, Issue 12, 
% December 2020, giaa133, https://doi.org/10.1093/gigascience/giaa133.
% The current script uses the code from the following github repository:
% https://github.com/nrdg/fracridge. We have included the toolbox in our external
% non-default package and use it in this script.
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
%     A matrix of the target variables to be predicted (dimension:
%     #subjects x K). K is the number of target variables to be predicted.
% 
%   - params.num_inner_folds (optional)
%     A string or a scalar. The number of inner-loop folds within each
%     training fold for hyperparameter selection.
%     If params does not have a field called num_inner_folds, this script
%     uses the same number of outer-loop folds as the number of inner-loop
%     folds.
% 
%   - params.lambda (optional)
%     A vector of the fraction parameters between 0 and 1. Fractions can be
%     exactly 0 or exactly 1. However, values in between 0 and 1 should be 
%     no less than 0.001 and no greater than 0.999.
%     For example, params.lambda could be 0:.05:1 or 0:.1:1. 
%     If this field is not passed in, the script will use the default lambda
%     values: params.lambda = [0.05:0.05:1]. 
%
%   - params.outdir (compulsory)
%     A string of the full-path output directory.
% 
%   - params.outstem (compulsory)
%     A string appended to the filename to specify the output files (y
%     after regression, accuracy files, ...). For example, if outstem =
%     'AngAffect_Unadj', then the final optimal accuracy file will be named 
%     as fullfile(outdir, 'results', 'optimal_acc', 'AngAffect_Unadj.mat'),
%     and the optimal model hyperparameters of k-th fold will be saved in
%     fullfile(outdir, 'params', ['dis_' k '_cv'], ...
%     'selected_parameters_AngAffect_Unadj.mat').
%
%   - params.metric (optional)
%     A string, specifying how the value of lambda will be optimised with.
%     Can choose the following: corr, COD, predictive_COD, MAE,
%     MAE_norm,MSE,MSE_norm. Default is predictive_COD
%
%  output (saved)
%   - acc_metric_train - training accuracy (given in correlation)
%   - acc_corr_test - test accuracies (given in correlation)
%   - y_predict - predicted target values
%   - optimal_statistics - a cell array of size equal to the number of
%   folds. Each element in the cell array is a structure storing the
%   accuracies of each possible accuracy metric (eg. corr, MAE, etc).
% 
% Written by Yanrui, Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Preparation
if(~isfield(params, 'num_inner_folds'))
    params.num_inner_folds = length(params.sub_fold);
end

if(~isfield(params, 'lambda'))
    params.lambda = [0.05:0.05:1];
end

if(~isfield(params, 'metric'))
    params.metric = 'predictive_COD';
end

%% step 1. Regress covariates from y
fprintf('# step 1: regress covariates from y for each fold.\n')
CBIG_crossvalid_regress_covariates_from_y( ...
    params.y, params.covariates, params.sub_fold, params.outdir, params.outstem);

%% step 2. For each fold, optimize the linear ridge regression model
param_dir = fullfile(params.outdir, 'params');
fprintf('# step 2: optimize linear ridge regression model and predict target measure.\n')
y_predict = cell(length(params.sub_fold), 1);
acc_corr_test = zeros(length(params.sub_fold), size(params.y,2));
for currfold = 1:length(params.sub_fold)
    curr_param_dir = fullfile(param_dir, ['fold_' num2str(currfold)]);
    if(~exist(curr_param_dir, 'dir'))
        mkdir(curr_param_dir)
    end
    fprintf('=========================== Fold %d ===========================\n', currfold)
    
    load(fullfile(params.outdir, 'y', ['fold_' num2str(currfold)], ['y_regress_' params.outstem '.mat']))
    y_train = y_resid(params.sub_fold(currfold).fold_index==0,:);
    y_test = y_resid(params.sub_fold(currfold).fold_index==1,:);
    
    if(ndims(params.feature_mat) == 3)
        params.feature_mat = reshape_3D_features(params.feature_mat);
    end

    feature_train = params.feature_mat(:, params.sub_fold(currfold).fold_index==0);
    feature_test = params.feature_mat(:, params.sub_fold(currfold).fold_index==1);
    
    %% step 2.1. select hyperparameters: inner-loop CV
    param_file = fullfile(curr_param_dir, ['selected_parameters_' params.outstem '.mat']);
    if(~exist(param_file, 'file'))
        [loss, y_pred, acc_corr_fold] = CBIG_LRR_frac_innerloop_cv( feature_train', y_train,...
         params.lambda, params.num_inner_folds, params.metric );
        
        loss_mean = mean(loss,1);
        [min_loss, min_idx] = min(loss_mean,[],3);
        curr_lambda = params.lambda(min_idx);
        acc_train_mean = mean(acc_corr_fold,1);
        acc_train_mean = reshape(acc_train_mean, size(acc_train_mean,2), size(acc_train_mean,3));
        curr_acc_train = acc_train_mean(sub2ind(size(acc_train_mean),1:size(acc_train_mean,1),min_idx));
        for i = 1:size(loss,1)
            curr_y_pred{i} = y_pred{i}(:,sub2ind(size(acc_train_mean),1:size(acc_train_mean,1),min_idx));
        end


        save(param_file, 'curr_lambda', 'min_loss')
        save(fullfile(curr_param_dir, ['acc_train_' params.outstem '.mat']), 'curr_acc_train', 'curr_y_pred');
    else

        load(param_file)
        if(~exist(fullfile(curr_param_dir, ['acc_train_' params.outstem '.mat']), 'file'))
            [~, curr_y_pred, curr_acc_train] = CBIG_LRR_frac_innerloop_cv( feature_train',...
             y_train, curr_lambda, params.num_inner_folds, params.metric );

            save(fullfile(curr_param_dir, ['acc_train_' params.outstem '.mat']), 'curr_acc_train', 'curr_y_pred');
        else
            load(fullfile(curr_param_dir, ['acc_train_' params.outstem '.mat']))
        end
    end
    opt_lambda(currfold, :) = curr_lambda;
    acc_metric_train(currfold, :) = curr_acc_train;
    disp('>>>>>>>>>>>>>>>>>>>>>>>>> Innerloop Done')

    %% Training & test
    fprintf('Test fold %d:\n', currfold);
    [acc_corr, optimal_statistics{currfold},y_predict{currfold}] = CBIG_LRR_frac_train_test( ...
    feature_train', feature_test', y_train, y_test, curr_lambda );
    
    acc_corr_test(currfold, :) = acc_corr;
    clear acc_corr
    disp('>>>>>>>>>>>>>>>>>>>>>>>>> Test Done')
end

opt_out = fullfile(params.outdir, 'results', 'optimal_acc', [params.outstem '.mat']);
mkdir(fullfile(params.outdir, 'results', 'optimal_acc'))
save(opt_out, 'acc_metric_train', 'acc_corr_test', 'y_predict', 'optimal_statistics')
fprintf('Finished!\n')

end


function out = reshape_3D_features(features)
    temp = ones(size(features,1), size(features,2));
    tril_ind = tril(temp, -1);
    
    features_reshaped = reshape(features, size(features,1)*size(features,2), size(features, 3));
    out = features_reshaped(tril_ind==1, :);
end

