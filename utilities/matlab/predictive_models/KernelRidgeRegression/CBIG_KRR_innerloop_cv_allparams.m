function [acc, loss] = CBIG_KRR_innerloop_cv_allparams( test_fold, sub_fold, num_inner_folds, ...
    data_dir, y_resid_stem, with_bias, ker_param, lambda_set, threshold_set, metric )

% [acc, acc_concat] = CBIG_KRR_innerloop_cv_allparams( test_fold, num_inner_folds, ...
%     outdir, y_resid, ker_param, lambda_set, threshold_set )
% 
% This function performs inner-loop cross-validation of kernel ridge
% regression with all hyperparameter combinations (differnt kernels,
% lambda, and/or threshold). The inner-loop cross-validation will be used
% to determine the optimal hyperparameter. This function runs inner-loop CV
% for a single test fold.
% 
% Note: if your target variables to be predicted, y, contain both binary
% (e.g. sex) and non-binary measures. You need to predict them separately
% (i.e. save the y values in different files with different "y_resid_stem",
% and run this function twice).
% 
% Inputs:
%   - test_fold
%     A string or a scalar, current test fold.
% 
%   - sub_fold
%     The data split for cross-validation.
%     It is a num_test_folds x 1 structure with a field "fold_index".
%     sub_fold(i).fold_index is a #subjects x 1 binary vector, where 1
%     refers to the test subjects in the i-th fold; 0 refers to the
%     training subjects in the i-th fold.
%     If the user does not need cross-validation, the length of sub_fold
%     will be 1. In this case, sub_fold.fold_index has 3 unique values: 
%     0 - training set;
%     1 - test set;
%     2 - validation set.
% 
%   - num_inner_fold
%     A string or a scalar, number of inner-loop folds. 
%     For cross-validation, each training set will be split into multiple
%     inner-loop folds for hyperparameter selection.
%     If the user has enough data and does not need cross-validation; and
%     the data were split into training, validation and test set, then
%     num_inner_fold = 1.
% 
%   - data_dir
%     Full path of the input/output directory.
%     A subfolder [data_dir '/innerloop_cv/fold_' num2str(test_fold)] will
%     be created to save the innerloop cross-validation results.
% 
%   - y_resid_stem
%     A string that was used to describe the .mat file of the target
%     variable, y, after regression (which was generated using
%     CBIG_crossvalid_regress_covariates_from_y.m). For example, if the
%     filename is <path_to_file>/y_regress_58behaviors.mat, "y_resid_stem"
%     will be '58behaviors'.
%     See the description of "outstem" parameter of function
%     CBIG_crossvalid_regress_covariates_from_y.m
% 
%   - with_bias
%     A a scalar (choose from 0 or 1).
%     - with_bias = 0 means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
%   
%   - ker_param (optional)
%     A K x 1 structure with two fields: type and scale. K denotes the
%     number of kernels.
%     ker_param(k).type contains a string that specifies the type of the k-th kernel. 
%                       Choose from
%                       'corr'        - Pearson's correlation;
%                       'Gaussian'    - Gaussian kernel;
%                       'Exponential' - exponential kernel.
%     ker_param(k).scale is a scalar specifying the scale of k-th kernel
%     (for Gaussian kernel or exponential kernel). If ker_param(k).type == 'corr',
%     ker_param(k).scale = NaN.
%     If "ker_param" is not passed in, only correlation kernel will be
%     calculated.
% 
%   - lambda_set (optional)
%     A vector of numbers for grid search of lambda (the regularization
%     parameter). If "lambda_set" is not passed in or is 'NONE', it will be
%     set as default:
%     [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
%        5 10 15 20]
% 
%   - threshold_set (optional)
%     A vector of numbers (thresholds). Different thresholds used will
%     affect the accuracy of the target variable (binary) predictions. To
%     be used when the target variable to predict is binary (eg. sex). The
%     threshold is used as a point to divide the prediction into the two
%     classes (eg. Male or Female). "Thresholds" should be between -1 and
%     1. 
%     If "threshold_set" is not passed in or is 'NONE', it will
%     be set as default: [-1:0.1:1].
%
%   - metric (optional)
%     A string indicating the metric used to define prediction loss. The
%     loss is used to choose hyperparameters.
%     Choose from:
%       'corr'              - Pearson's correlation;
%       'COD'               - Coefficient of determination. Defined as
%                             1-||y_pred-y_test||^2/||mean(y_test)-y_test||^2,
%                             where y_pred is the prediction of the test data, 
%                             y_test is the groud truth of the test data, 
%                             and mean(y_test) is the mean of test data
%       'predictive_COD'    - Predictive coefficient of determination. Defined as
%                             1-||y_pred-y_test||^2/||mean(y_train)-y_test||^2,
%                             where y_pred is the prediction of the test data, 
%                             y_test is the groud truth of the test data, 
%                             and mean(y_train) is the mean of training data
%       'MAE'               - mean absolute error
%       'MAE_norm'          - mean absolute error divided by the standard
%                             derivation of the target variable of the training set
%       'MSE'               - mean squared error
%       'MSE_norm'          - mean squared error divided by the variance
%                             of the target variable of the training set
% 
% 
% Outputs:
%   - acc
%     A cell with dimension of #kernels x #lambda x #thresholds of
%     inner-loop CV accuracy for each hyperparameter combination, averaged
%     across inner-loop folds. Each cell array is a 1 x #TargetVariables
%     vector.
%  
%   - loss
%     A cell with dimension of #kernels x #lambda x #thresholds of
%     inner-loop CV loss for each hyperparameter combination, averaged
%     across inner-loop folds. Each cell array is a 1 x #TargetVariables
%     vector.
%    
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Jingwei Li and Ru(by) Kong

%% setting up
if(ischar(test_fold))
    test_fold = str2num(test_fold);
end

if(ischar(num_inner_folds))
    num_inner_folds = str2num(num_inner_folds);
end

if(~isempty(y_resid_stem))
    y_resid_stem = ['_' y_resid_stem];
end
y_resid_file = fullfile(data_dir, 'y', ['fold_' num2str(test_fold)], ...
    ['y_regress' y_resid_stem '.mat']);
load(y_resid_file)

% set default hyperparamters if not passed in
if(~exist('ker_param', 'var') || strcmpi(ker_param, 'none'))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
end

if(~exist('lambda_set', 'var') || strcmpi(lambda_set, 'none'))
    lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
        5 10 15 20];
end

if(~exist('metric', 'var') || strcmpi(metric, 'none'))
    metric = 'predictive_COD';
end

for i = 1:size(y_orig, 2)
    if(length(unique(y_orig(:,i))) > 2)
        % not binary
        bin_flag(i) = 0;
    else
        % binary case
        bin_flag(i) = 1;
    end
end
if(any(bin_flag==1) && any(bin_flag==0))
    error('Mixture of binary (e.g. sex) and continuous cases. Please run them separately.')
elseif(~any(bin_flag==1))    % all continuous
    bin_flag = 0;
    threshold_set = NaN;
else                         % all binary
    bin_flag = 1;
    if(~exist('threshold_set', 'var') || strcmpi(threshold_set, 'none') || isempty(threshold_set))
        threshold_set = [-1:0.1:1];
    end
end

% check if using cross-validation
if(length(sub_fold)>1)
    % cross-validation
    num_valid_sub = [];
    kernel_dir = [data_dir '/FSM_innerloop/fold_' num2str(test_fold)];
    sub_index = sub_fold(test_fold).fold_index==0;
else
    num_valid_sub = length(find(sub_fold.fold_index==2));
    kernel_dir = [data_dir '/FSM_innerloop'];
    % first S1 subjects are from training set, 
    % next S2 subject are from validation set
    sub_index = [find(sub_fold.fold_index==0); find(sub_fold.fold_index==2)];
end

%% inner-loop cv for each parameter combination
outdir = fullfile(data_dir, 'innerloop_cv', ['fold_' num2str(test_fold)]);

if(~exist(fullfile(outdir, ['acc' y_resid_stem '.mat']), 'file'))
    curr_y_resid = y_resid(sub_index, :);
    curr_y_orig = y_orig(sub_index, :);
    
    for k = 1:length(ker_param)
        fprintf('Kernel type: %s, scale: %f\n', ker_param(k).type, ker_param(k).scale)
        if(strcmp(ker_param(k).type, 'corr'))
            kernel = fullfile(kernel_dir, ['FSM_' ker_param(k).type '.mat']);
        else
            kernel = fullfile(kernel_dir, ['FSM_' ker_param(k).type '_' num2str(ker_param(k).scale) '.mat']);
        end
        load(kernel)
        
        for l = 1:length(lambda_set)
            lambda = lambda_set(l);
            fprintf('  lambda: %f\n', lambda);
            
            for t = 1:length(threshold_set)
                threshold = threshold_set(t);
                fprintf('    threshold: %f\n', threshold)
                [y_p, y_t, acc_fold, loss_fold] = CBIG_KRR_innerloop_cv( bin_flag, num_inner_folds, ...
                    FSM, curr_y_resid, curr_y_orig, num_valid_sub, with_bias, lambda, threshold, metric );
                clear threshold
                y_pred{k,l,t} = y_p;
                acc{k,l,t} = acc_fold;
                loss{k,l,t} = loss_fold;
            end
            clear lambda
        end
        clear FSM
    end
    
    if(~exist(outdir, 'dir'))
        mkdir(outdir)
    end
    save(fullfile(outdir, ['acc' y_resid_stem '.mat']), 'acc', 'loss', 'y_pred');
else
    load(fullfile(outdir, ['acc' y_resid_stem '.mat']))
    if(size(acc,1)~=length(ker_param) || size(acc,2)~=length(lambda_set) || size(acc,3)~=length(threshold_set))
        error('Wrong size of ''acc'' in existing file %s. Please remove it and rerun this script.', ...
            [outdir '/acc' y_resid_stem '.mat'])
    end
end

end

