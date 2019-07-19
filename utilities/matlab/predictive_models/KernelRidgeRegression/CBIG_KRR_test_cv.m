function [y_p, y_t, acc] = CBIG_KRR_test_cv( bin_flag, kernel_train, kernel_test, ...
    y_resid_train, y_resid_test, y_orig_test, with_bias, lambda, threshold )

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
% Written by Jingwei Li, Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% setting up
if(bin_flag==1 && (~exist('y_orig_test', 'var') || isempty(y_orig_test)) )
    error('"y_orig_test" variable is needed when y is binary (e.g. sex)')
end

%%
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
        X = ones(N,1);
        beta{i} = (X' * (Kr\ X)) \ X' * (Kr \ y);
        alpha{i} = Kr \ (y - X * beta{i});
    end
    
    %% Test
    K_test = kernel_test(~isnan(y_resid_test(:,i)), ~nan_index);
    y_p{i} = nan(size(y_resid_test, 1), 1);
    if(with_bias==0)
        %%%%%%%%%%%%%%%%%%%% Without bias term
        y_out = K_test * alpha{i};
        y_p{i}(~isnan(y_resid_test(:,i))) = y_out;
    else
        %%%%%%%%%%%%%%%%%%%% With bias term
        N_test = sum(~isnan(y_resid_test(:,i)));
        y_out = K_test * alpha{i} + ones(N_test,1) .* beta{i};
        y_p{i}(~isnan(y_resid_test(:,i))) = y_out;
    end
    
    if(bin_flag==1)
        y_t{i} = y_orig_test(~isnan(y_orig_test(:,i)), i);
        TP = length(find((y_out > threshold) & (y_t{i}==1) == 1));
        TN = length(find((y_out < threshold) & (y_t{i}==0) == 1));
        acc(i) = (TP + TN) / length(y_t{i});
        y_t{i} = y_orig_test(:,i);
    else
        y_t{i} = y_resid_test(~isnan(y_resid_test(:,i)), i);
        acc(i) = CBIG_corr(y_out, y_t{i});
        y_t{i} = y_resid_test(:, i);
    end
end


end

