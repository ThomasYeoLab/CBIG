function [acc, optimal_stats,pred] = ...
    CBIG_LRR_frac_train_test( feat_train, feat_test, y_train, y_test, lambda)

% [acc, optimal_stats,pred,pred_train] = ...
%    CBIG_LRR_frac_train_test( feat_train, feat_test, y_train, y_test, lambda, metric )
% 
% Based on the specified fraction hyperparameter (lambda), this function 
% trains a fractional linear ridge regression model to predict the target
% value (y_train) from the features (feat_train) in the training set. After
% that, it applies the trained model on the test set (y_test and feat_test).
% This function is an adaptation from CBIG_LRR_train_test in CBIG_repo. It
% uses the fitridge function provided by: https://github.com/nrdg/fracridge.
% 
% Inputs:
%   - feat_train:
%     A 2D matrix (#training subjects x #features). It is the selected
%     features of training subjects.
% 
%   - feat_test:
%     A 2D matrix (#test subject x x#features). It is the selected features
%     of test subjects.
%     
%   - y_train:
%     A matrix of the target measures in the training set to be predicted
%     (dimension: #subjects x K). K is the number of target variables to be
%     predicted.
% 
%   - y_test:
%     A matrix of the target measures in the test set to be predicted
%     (dimension: #subjects x K). K is the number of target variables to be
%     predicted.
% 
%   - lambda:
%     A scalar, the fraction hyperparameter. Value is between 0 and 1.
%
% Outputs:
%   - acc:
%     A vector of the prediction accuracy of the test set. The length of the
%     vector is the number of target variables to be predicted (dimension:
%     1 x K). K is the number of target variables to be predicted.
%
%   - optimal_stats:
%     A structure containing the prediction accuracy of the test set
%     using different accuracy metric (eg. corr, MAE, etc).
%
%   - pred:
%     A matrix of the predicted target measures in the test set (dimension:
%     #subjects x K). K is the number of target variables to be predicted.
% 
% Written by Yanrui, Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
        'non_default_packages', 'fracridge'));

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

p = size(feat_train, 2);

%% training

beta_tmp = fracridge(feat_train,lambda,y_train); %beta_tmp is a #feature x #lambda x #behaviors
for b = 1:size(beta_tmp,3)
    beta(:,b) = beta_tmp(:,b,b);
end

%% Predicting
pred = feat_test * beta; %pred: #sub x #behaviors

for i = 1:length(metrics)
    curr_metric = metrics{i};
    for b = 1:size(pred,2)
        [acc(b), ~] = CBIG_compute_prediction_acc_and_loss(pred(:,b),y_test(:,b), curr_metric,y_train(:,b));
    end
    optimal_stats.(curr_metric) = acc;
end

acc = optimal_stats.corr;

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
'non_default_packages', 'fracridge'));

end