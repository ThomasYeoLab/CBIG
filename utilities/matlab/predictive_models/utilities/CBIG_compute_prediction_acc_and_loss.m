function [pred_stats,loss] = CBIG_compute_prediction_acc_and_loss(y_pred,y_test,metric,y_train)

% [pred_stats,loss] = CBIG_compute_prediction_acc_and_loss(y_pred,y_test,metric,y_train)
%
% This function computes the prediction accuracy and loss defined by the
% metric passed in.
%
% Inputs:
%   - y_pred
%     Vector of length # test subjects. Prediction of the traget variable of
%     the test subjects
%
%   - y_test
%     Vector of length # test subjects. Groudtruth of the target variable of
%     the test subjects
%
%   - metric
%     A string indicating the prediction stats to be computed.
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
%   - y_train
%     Vector of length # test subjects. The target variable of the traning 
%     subjects. Only useful when compute the predictive_COD, MAE_norm, or MSE_norm

% Outputs:
%   - pred_stats
%     A scalar. Prediction statistics of the given metric
%
%   - loss
%     A scalar. Prediction loss of the given the metric. The loss is the same 
%     as pred_stats if the given metric is the smaller the better (e.g. 'MAE','MSE').
%     The loss is the addictive inverse of the pred_stats if the given metric is
%     the bigger the better (e.g. 'corr', 'COD').
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

switch metric
    case 'corr'
        pred_stats = CBIG_corr(y_test,y_pred);
        loss = -pred_stats;
    case 'COD'
        ss_res = sum((y_pred - y_test).^2);
        ss_total = sum((y_test - mean(y_test)).^2);
        pred_stats = 1-ss_res/ss_total;
        loss = -pred_stats;
    case 'predictive_COD'
        ss_res = sum((y_pred - y_test).^2);
        ss_total = sum((y_test - mean(y_train)).^2);
        pred_stats = 1-ss_res/ss_total;
        loss = -pred_stats;
    case 'MAE'
        pred_stats = mean(abs(y_pred-y_test));
        loss = pred_stats;
    case 'MAE_norm'
        pred_stats = mean(abs(y_pred-y_test))/std(y_train);
        loss = pred_stats;
    case 'MSE'
        pred_stats = mean((y_pred - y_test).^2);
        loss = pred_stats;
    case 'MSE_norm'
        pred_stats = mean((y_pred - y_test).^2)/var(y_train);
        loss = pred_stats;
    otherwise
        disp(metric);
        error('Unexpected metric');
end