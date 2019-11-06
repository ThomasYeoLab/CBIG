function [acc_corr, mse] = CBIG_LRR_innerloop_cv( feature_mat, y, lambda, tot_folds )

% [acc_corr, mse] = CBIG_LRR_innerloop_cv( feature, y, lambda, tot_folds )
% 
% This function performs the inner-loop cross-validation of linear ridge
% regression for a given lambda and selected feature set.
% 
% Inputs:
%   - feature_mat:
%     A 2D matrix of the selected features (dim: #subjects x #features).
%     This function automatically prepends the bias term (a vector of ones)
%     onto "feature_mat".
% 
%   - y:
%     A column vector of target measure to predict (dim: #subjects x 1).
% 
%   - lambda:
%     A scalar, the regularization hyperparameter.
% 
%   - tot_folds:
%     A scalar, how many inner-loop CV folds the user wants to conduct.
% 
% Outputs:
%   - acc_corr:
%     A scalar. It is the mean correlation accuracy averaged across
%     inner-loop folds.
% 
%   - mse:
%     A scalar. It is the mean mean sqaured error averaged across
%     inner-loop folds.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(size(y,2)~=1)
    error('''y''should be a column vector.')
end


%% split data
rng('default');  rng(1);
num_sub = size(feature_mat, 1);
feature_mat = cat(2, ones(num_sub, 1), feature_mat);
p = size(feature_mat, 2);
cv_idx = cvpartition(num_sub, 'kfold', tot_folds);


%% For each fold, training and testing
for curr_fold = 1:tot_folds
    
    curr_train_idx = cv_idx.training(curr_fold);
    curr_test_idx = cv_idx.test(curr_fold);
    
    feat_train = feature_mat(curr_train_idx, :);
    feat_test = feature_mat(curr_test_idx, :);
    
    y_train = y(curr_train_idx);
    y_test = y(curr_test_idx);
    
    % training
    beta = (feat_train' * feat_train + lambda*eye(p))' \ feat_train' * y_train;
    y_predict_train = feat_train * beta;
    acc_corr_train = CBIG_corr(y_train, y_predict_train);
    mse_train = mean((y_train - y_predict_train).^2);
    
    % test
    y_predict = feat_test * beta;
    acc_corr_fold(curr_fold, 1) = CBIG_corr(y_test, y_predict);
    mse_fold(curr_fold, 1) = mean((y_test - y_predict).^2);
end

acc_corr = mean(acc_corr_fold);
mse = mean(mse_fold);
fprintf('Mean testing accuracy (correlation) across %d folds: %f; mean mse: %f\n', tot_folds, acc_corr, mse);

end

