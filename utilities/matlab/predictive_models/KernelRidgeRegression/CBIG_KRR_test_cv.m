function [y_p, y_t, acc, pred_stats, y_pred_train] = CBIG_KRR_test_cv( bin_flag, kernel_train, kernel_test, ...
    y_resid_train, y_resid_test, y_orig_test, with_bias, lambda, threshold, saving_stats )

% [y_p, y_t, acc] = CBIG_KRR_test_cv( bin_flag, kernel_train, kernel_test, ...
%     y_resid_train, y_resid_test, y_orig_test, lambda, threshold )
%
% This function computes the test accuracies of all target variables (y) in
% a test fold given a particular hyperparameter combination of kernel
% parameter, regularization parameter, and/or threshold.
%
% Inputs:
%   - bin_flag
%     0 or 1. 1 means the target variables to be predicted (i.e. y_orig)
%     are binary. 0 means none of them are binary.
%     If the target variables consist of both binary and non-binary
%     measures, you need to predict them separately with different
%     bin_flag.
%
%   - kernel_train
%     A #TrainingSubject x #TrainingSubject kernel matrix. The ordering of
%     subjects should match the ordering of y_resid_train.
%     It is supposed to be an output of CBIG_KRR_generate_kernels.m.
%
%   - kernel_test
%     A #TestSubject x #TrainingSubject kernel matrix. The ordering of
%     subjects should match the ordering of y_resid_test.
%     It is supposed to be an output of CBIG_KRR_generate_kernels.m
%
%   - y_resid_train
%     A #TrainingSubject x #TargetVariable matrix of y after nuisance
%     regression. It is supposed to be an output of
%     CBIG_crossvalid_regress_covariates_from_y.m
%
%   - y_resid_test
%     A #TestSubject x #TargetVariable matrix of y after nuisance
%     regression. It is supposed to be an output of
%     CBIG_crossvalid_regress_covariates_from_y.m
%
%   - y_orig_test (optional)
%     A #TestSubject x #TargetVariable matrix of y before nuisance
%     regression. This input is required for the calculation of prediction
%     accuracy when the target variable y is binary (e.g. sex).
%     If y is not a binary variable to predict, specify y_orig = [].
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
%     A scalar (threshold). Thresholds are used when the target variable to
%     predict is binary (eg. sex). The threshold is used as a point
%     (between -1 and 1) to divide the prediction into the two classes (eg.
%     Male or Female).
%     Ignore this parameter if y is not binary.
%   - saving_stats
%     A cell array of strings. Each element of the cell array indicates the
%     prediction statistics you want to compute and save.
%     Supported metric:
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
%     This argument should be a cell array of one or multiple strings
%     listed above. 
%     Default:{'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
%
% Outputs:
%   - y_p
%     Cell arrays (length #TargetVariable) of predicted target variables y.
%     Each cell array is an #TestSubject x 1 vector storing the predicted
%     values of all the test subjects.
%
%   - y_t
%     Cell arrays (length #TargetVariable) of groud truth of target
%     variables y. Each cell array is an #TestSubject x 1 vector storing the
%     ground truth values of all the test subjects.
%     If y is binary (e.g. sex), y_t will be the original y before
%     regression of nuisace covariates. If y is not binary, y_t will be the
%     y after regression of nuisance covariates.
%
%   - acc
%     A 1 x size(y_resid,2) vector, the accuracy of each target variable
%     prediction.
%
%   - pred_stats
%     A length(saving_stats) X #TargetVariable matrix. The prediction statistics
%     for each metric defined in the input argument saving_stats for each
%     target variable.
%
%   - y_pred_train
%     Cell arrays (length #TargetVariable) of predicted target variables y of the 
%     training subjects. It can be used for model interpretation.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Jingwei Li and Ru(by) Kong

%% setting up
if(bin_flag==1 && (~exist('y_orig_test', 'var') || isempty(y_orig_test)) )
    error('"y_orig_test" variable is needed when y is binary (e.g. sex)')
end

%%
pred_stats = zeros(length(saving_stats),size(y_resid_train, 2));
if sum(sum(isnan(y_resid_train))) > 0 || sum(sum(isnan(y_resid_test))) >0
    for i = 1:size(y_resid_train, 2)
        %% Training
        nan_index = isnan(y_resid_train(:,i));
        K = kernel_train(~nan_index,~nan_index);
        N = sum(~nan_index);
        y = y_resid_train(~nan_index,i);
        Kr = K + lambda * eye(N);
        if(with_bias==0)
            %%%%%%%%%%%%%%%%%%%% Without bias term
            % alpha = (K + lambda * I)^-1 * y
            %  Nx1    NxN   1x1    NxN      Nx1
            alpha{i} = Kr \ y;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            %%%%%%%%%%%%%%%%%%%% With bias term
            inv_Kr = inv(Kr);
            X = ones(N,1);
            beta{i} = (X' * (inv_Kr * X)) \ X' * (inv_Kr * y);
            alpha{i} = inv_Kr * (y - X * beta{i});
        end
        
        %% Test
        K_test = kernel_test(:, ~nan_index);
        if(with_bias==0)
            %%%%%%%%%%%%%%%%%%%% Without bias term
            y_out = K_test * alpha{i};
            y_p{i} = y_out;
            y_pred_train{i} = NaN(length(nan_index),1);
            y_pred_train{i}(~nan_index) = K * alpha{i};
        else
            %%%%%%%%%%%%%%%%%%%% With bias term
            N_test = size(y_resid_test, 1);
            y_out = K_test * alpha{i} + ones(N_test,1) .* beta{i};
            y_p{i} = y_out;
            y_pred_train{i} = NaN(length(nan_index),1);
            y_pred_train{i}(~nan_index) = K * alpha{i} + ones(N,1) .* beta{i};
        end
        y_out = y_out(~isnan(y_resid_test(:,i)));
        
        if(bin_flag==1)
            y_t{i} = y_orig_test(~isnan(y_orig_test(:,i)), i);
            TP = length(find((y_out > threshold) & (y_t{i}==1) == 1));
            TN = length(find((y_out < threshold) & (y_t{i}==0) == 1));
            acc(i) = (TP + TN) / length(y_t{i});
            y_t{i} = y_orig_test(:,i);
            pred_stats(:,i) = nan;
        else
            y_t{i} = y_resid_test(~isnan(y_resid_test(:,i)), i);
            acc(i) = CBIG_corr(y_out, y_t{i});
            for metric_ind = 1:length(saving_stats)
                [pred_stats(metric_ind,i),~] = CBIG_compute_prediction_acc_and_loss(...
                    y_out, y_t{i}, saving_stats{metric_ind}, y);
            end
            y_t{i} = y_resid_test(:, i);            
        end
    end
else
    K = kernel_train;
    N = size(K,1);
    y = y_resid_train;
    Kr = K + lambda * eye(N);
    if(with_bias==0)
        %%%%%%%%%%%%%%%%%%%% Without bias term
        % alpha = (K + lambda * I)^-1 * y
        %  NxV    NxN   1x1    NxN     NxV (V:# target variables)
        alpha = Kr \ y;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    else
        %%%%%%%%%%%%%%%%%%%% With bias term
        inv_Kr =  inv(Kr);
        X = ones(N,1);
        beta = (X' * (inv_Kr * X)) \ X' * (inv_Kr * y);
        alpha = inv_Kr * (y - X * beta);
    end
    
    %% Test
    K_test = kernel_test;
    if(with_bias==0)
        %%%%%%%%%%%%%%%%%%%% Without bias term
        y_out = K_test * alpha;
        y_train_out = K * alpha;
    else
        %%%%%%%%%%%%%%%%%%%% With bias term
        N_test = size(y_resid_test,1);
        y_out = K_test * alpha + ones(N_test,1) .* beta;
        y_train_out = K * alpha + ones(N,1) .* beta;
    end
    
    for i = 1:size(y_resid_train,2)
        y_p{i} = y_out(:,i);
        y_pred_train{i} = y_train_out(:,i);
        if(bin_flag==1)
            y_t{i} = y_orig_test(:,i);
            TP = length(find((y_out(:,i) > threshold) & (y_t{i}==1) == 1));
            TN = length(find((y_out(:,i) < threshold) & (y_t{i}==0) == 1));
            acc(i) = (TP + TN) / length(y_t{i});
            pred_stats(:,i) = nan;
        else
            y_t{i} = y_resid_test(:, i);
            acc(i) = CBIG_corr(y_out(:,i), y_t{i});
            for metric_ind = 1:length(saving_stats)
                [pred_stats(metric_ind,i),~] = CBIG_compute_prediction_acc_and_loss(...
                    y_out(:,i), y_t{i}, saving_stats{metric_ind}, y(:,i));
            end
        end
    end
end


end