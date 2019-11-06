function [y_p, y_t, acc, loss] = CBIG_KRR_innerloop_cv( bin_flag, num_inner_folds, ...
    kernel, y_resid, y_orig, num_valid_sub, with_bias, lambda, threshold, metric )

% [y_p, y_t, acc, acc_concat] = CBIG_KRR_innerloop_cv( bin_flag, num_inner_folds, ...
%     kernel, y_resid, y_orig, num_valid_sub, lambda, threshold )
%
% This function performs innner-loop cross-validation of kernel ridge
% regression with only one hyperparamter combination (kernel type and
% scale, regularization, threshold for binary cases (if necessary, e.g.
% sex)).
%
% Inputs:
%   - bin_flag
%     0 or 1. 1 means the target variables to be predicted (i.e. y_orig)
%     are binary. 0 means none of them are binary.
%     If the target variables consist of both binary and non-binary
%     variables, you need to call this function twice and predict them
%     separately with different "bin_flag".
%
%   - num_inner_folds
%     A scalar, the number of inner-loop cross-validation folds for
%     hyperparameter selection.
%     If there are enough data so that cross-validation is not needed, set
%     num_inner_folds = 1.
%
%   - kernel
%     An S x S kernel matrix.
%     For cross-validation case, S is the number of training subjects. The
%     ordering of training subjects needs to be matched with the ordering
%     of the target variable y.
%     If data were split into training, validation, and test sets, then S
%     is the number of training subject + the number of validation
%     subjects. The first S1 subjects of S will be the training subjects
%     and the last S2 = S - S1 subjects will be the validation subjects
%
%   - y_resid
%     An S x #VariableToPredict matrix of the target variable y after
%     regressing out nuisance regressors.
%     If cross-validation is not used, the first S1 rows of y_resid will
%     correspond to the training subjects and the last S2 = S - S1 rows of
%     y_resid will correspond to the validation subjects.
%
%   - y_orig (optional)
%     An S x #VariableToPredict vector of the target variable y before
%     regression of covariates. This vector is needed for binary (0 or 1) case (e.g. sex)
%     for accuracy calculation.
%     If y is not a binary variable, specify y_orig = [].
%     If cross-validation is not used, the first S1 rows of y_orig will
%     correspond to the training subjects and the last S2 = S - S1 rows of
%     y_orig will correspond to the validation subjects.
%
%   - num_valid_sub (optional)
%     A scalar, the number of validation subjects for non-cross-validation
%     case. It is assumed that the last "num_valid_sub" subjects in
%     "kernel", "y_resid" and "y_orig" are the validation subjects.
%     For cross-validation case, this argument is not needed, set
%     num_valid_sub = [].
%
%   - with_bias
%     A scalar (choose from 0 or 1).
%     - with_bias = 0 means the algorithm is to minimize
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
%
%   - lambda
%     A scalar, the l2-regularization parameter.
%
%   - threshold (optional)
%     A scalar, it is a threshold required to determine the prediction (1
%     or 0) for binary target variables (e.g. sex). Ignore this parameter
%     if y is not binary.
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
% Outputs:
%   - y_p
%     A cell array of length #TargetVariables. Each cell array contains the
%     predicted target values of all training subjects (using the
%     inner-loop CV manner) or validation subjects (using the
%     training-validation-test manner). Within each cell, the predicted
%     values follow the same ordering of y_resid.
%     For cross-validation case, each cell array is an S x 1 vector.
%     For training-validation-test case, each cell array is an S2 x 1 vector.
%
%   - y_t
%     A cell array of length #TargetVariables. Each cell array contains the
%     groud truth (actual) target values of all training subjects (using
%     the inner-loop CV manner) or validation subject (using
%     training-validation-test manner). Within each cell, the ground-truth
%     values follow the same ordering of y_resid.
%     If y is binary (e.g. sex), y_t will be the original y before
%     regression of nuisace covariates. If y is not binary, y_t will be the
%     y after regression of nuisance covariates.
%     For cross-validation case, each cell array is an S x size(y_resid, 2)
%     matrix.
%     For training-valdation-test case, each cell array is an S2 x
%     size(y_resid, 2) matrix.
%
%   - acc
%     A 1 x size(y_resid,2) vector, the accuracy of the prediction of each
%     target variable, averaged across all inner-loop folds.
%     For binary case, acc are the average of (true positive & negative) /
%     #subjects; for continuous variables, acc measures how correlated the
%     predicted and actual values are.
%
%   - acc_concat
%     A 1 x size(y_resid,2) vector, the accuracy of each variable. the
%     accuracy is computed based on the concatenated predicted scores of
%     all inner-loop folds.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Jingwei Li and Ru(by) Kong

%% setting up
num_sub = size(kernel, 1);

if(bin_flag==1 && (~exist('y_orig', 'var') || isempty(y_orig)) )
    error('"y_orig" variable is needed when y is binary (e.g. sex)')
end

% training and test indices for each inner-loop fold
if(num_inner_folds == 1)
    % all subjects were split into training, validation and test sets
    assert(~isempty(num_valid_sub), ...
        '''num_valid_sub'' needs to be passed in since cross-validation is not performed.')
    fold_index.training = 1:(num_sub - num_valid_sub);
    fold_index.test = (num_sub - num_valid_sub + 1):num_sub;
else
    % cross-validation case
    rng(1, 'twister');
    cv_idx = cvpartition(num_sub, 'kfold', num_inner_folds);
    
    for fold = 1:num_inner_folds
        fold_index(fold).training = cv_idx.training(fold);
        fold_index(fold).test = cv_idx.test(fold);
    end
end

%% kernel regression for each inner-loop fold
acc_fold = cell(num_inner_folds,size(y_resid, 2));
loss_fold = cell(num_inner_folds,size(y_resid, 2));
for fold = 1:num_inner_folds
    curr_train_idx = fold_index(fold).training;
    curr_test_idx = fold_index(fold).test;
    
    curr_ker_train = kernel(curr_train_idx, curr_train_idx);
    curr_ker_test = kernel(curr_test_idx, curr_train_idx);
    
    curr_y_train = y_resid(curr_train_idx, :);
    curr_y_test = y_resid(curr_test_idx, :);
    if(bin_flag==1)
        curr_y_true = y_orig(curr_test_idx, :);
    end
    
    if sum(sum(isnan(y_resid))) > 0
        for i = 1:size(y_resid, 2)
            %% Training
            nan_index = isnan(curr_y_train(:,i));
            K = curr_ker_train(~nan_index,~nan_index);
            N = sum(~nan_index);
            y = curr_y_train(~nan_index,i);
            Kr = K + lambda * eye(N);
            if(with_bias==0)
                %%%%%%%%%%%%%%%%%%%% Without bias term
                % alpha = (K + lambda * I)^-1 * y
                %  Nx1    NxN   1x1    NxN      Nx1
                alpha{fold,i} = Kr \ y;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            else
                %%%%%%%%%%%%%%%%%%%% With bias term
                inv_Kr = inv(Kr);
                X = ones(N,1);
                beta{fold,i} = (X' * (inv_Kr * X)) \ X' * (inv_Kr * y);
                alpha{fold,i} = inv_Kr * (y - X * beta{fold,i});
            end
            
            
            %% Test
            K_val = curr_ker_test(~isnan(curr_y_test(:,i)), ~nan_index);
            if(with_bias==0)
                %%%%%%%%%%%%%% Without bias term
                y_predict{fold,i} = K_val * alpha{fold,i};
            else
                %%%%%%%%%%%%%% With bias term
                N_val = sum(~isnan(curr_y_test(:,i)));
                y_predict{fold,i} = K_val * alpha{fold,i} + ones(N_val,1) .* beta{fold,i};
            end
            
            if(bin_flag==1)
                y_true{fold, i} = curr_y_true(~isnan(curr_y_true(:,i)), i);
                TP = length(find((y_predict{fold,i} > threshold) & (y_true{fold,i}==1) == 1));
                TN = length(find((y_predict{fold,i} < threshold) & (y_true{fold,i}==0) == 1));
                acc_fold{fold,i} = (TP + TN) / length(y_true{fold,i});
                loss_fold{fold,i} = -acc_fold{fold,i};
            else
                y_true{fold, i} = curr_y_test(~isnan(curr_y_test(:,i)), i);
                acc_fold{fold,i} = CBIG_corr(y_predict{fold,i}, y_true{fold,i});
                [~,loss_fold{fold,i}] = CBIG_compute_prediction_acc_and_loss(y_predict{fold,i},y_true{fold,i},metric,y);
            end
        end
    else
        %% Training
        K = curr_ker_train;
        y = curr_y_train;
        N = size(K,1);
        Kr = K + lambda * eye(N);
        if(with_bias==0)
            %%%%%%%%%%%%%%%%%%%% Without bias term
            % alpha = (K + lambda * I)^-1 * y
            %  NxV    NxN   1x1    NxN     NxV (V:# target variables)
            alpha = Kr \ y;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            %%%%%%%%%%%%%%%%%%%% With bias term
            inv_Kr = inv(Kr);
            X = ones(N,1);
            beta = (X' * (inv_Kr * X)) \ X' * (inv_Kr * y);
            alpha = inv_Kr * (y - X * beta);
        end
        
        
        %% Test
        K_val = curr_ker_test;
        if(with_bias==0)
            %%%%%%%%%%%%%% Without bias term
            y_predict_curr_fold = K_val * alpha;
        else
            %%%%%%%%%%%%%% With bias term
            N_val = size(curr_y_test,1);
            y_predict_curr_fold = K_val * alpha + ones(N_val,1) .* beta;
        end
        %{
        disp(size(y_predict_curr_fold));
        disp(size(K_val));
        disp(size(beta));
        %}
        for i = 1:size(y_resid, 2)
            y_predict{fold,i} = y_predict_curr_fold(:,i);
            if(bin_flag==1)
                y_true{fold, i} = curr_y_true(:,i);
                TP = length(find((y_predict{fold,i} > threshold) & (y_true{fold,i}==1) == 1));
                TN = length(find((y_predict{fold,i} < threshold) & (y_true{fold,i}==0) == 1));
                acc_fold{fold,i} = (TP + TN) / length(y_true{fold,i});
                loss_fold{fold,i} = -acc_fold{fold,i};
            else
                y_true{fold, i} = curr_y_test(:,i);
                acc_fold{fold,i} = CBIG_corr(y_predict{fold,i}, y_true{fold,i});
                [~,loss_fold{fold,i}] = CBIG_compute_prediction_acc_and_loss(y_predict{fold,i}...
                    ,y_true{fold,i},metric,y(:,i));
            end
        end      
    end
end

%% concatenate
for i = 1:size(y_resid, 2)
    y_predict_concat = [];
    y_true_concat = [];
    acc_all = 0;
    loss_all = 0;
    for fold = 1:num_inner_folds
        y_predict_concat = [y_predict_concat; y_predict{fold,i}];
        y_true_concat = [y_true_concat; y_true{fold,i}];
        acc_all = acc_all +  acc_fold{fold,i};
        loss_all = loss_all + loss_fold{fold,i};
    end
    acc(i) = acc_all/num_inner_folds;
    loss(i) = loss_all/num_inner_folds;
    y_p{i} = y_predict_concat;
    y_t{i} = y_true_concat;
end


end