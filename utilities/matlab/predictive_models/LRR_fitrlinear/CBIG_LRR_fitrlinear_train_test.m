function [acc, optimal_stats,pred,pred_train] = ...
    CBIG_LRR_fitrlinear_train_test( feat_train, feat_test, y_train, y_test, lambda, metric )

% [acc, optimal_stats,pred,pred_train] = ...
%    CBIG_LRR_fitrlinear_train_test( feat_train, feat_test, y_train, y_test, lambda, metric )
% 
% Based on the specified regularization hyperparameter (lambda), this
% function trains a linear ridge regression model to predict the target
% value (y_train) from the features (feat_train) in the training set. After
% that, it applies the trained model on the test set (y_test and feat_test).
% This function is an adaptation from CBIG_LRR_train_test in CBIG_repo. It
% uses the MATLAB function `fitrlinear` for fitting the LRR. The fitrlinear
% function automatically fits a bias term from input features. Hence, there
% is no need for the user to concatenate a bias term manually.
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
%     A column vector of the target measure in the training set (dim:
%     #training subjects x 1).
% 
%   - y_test:
%     A column vector of the target measure in the test set (dim: #test
%     subjects x 1).
% 
%   - lambda:
%     A scalar, the regularization hyperparameter.
%
%   - metric
%     A string, what is the accuracy metric that we want to optimise our
%     hyperparameter with. Can choose from
%     'corr','COD','predictive_COD','MAE_norm','MAE','MSE','MSE_norm'
% 
% Written by Yanrui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(size(y_train,2) > 1)
    error('''y_train'' should be a column vector')
end
if(size(y_test,2) > 1)
    error('''y_test'' should be a column vector')
end

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

p = size(feat_train, 2);
rng('default')
rng(1)

%% training
Mdl = fitrlinear(feat_train',y_train,'ObservationsIn','columns', 'Lambda',lambda, 'Learner',...
    'leastsquares', 'Regularization','ridge');

%% Predicting
pred = predict(Mdl,feat_test','ObservationsIn','columns');
pred_train = predict(Mdl,feat_train','ObservationsIn','columns');

for i = 1:length(metrics)
    curr_metric = metrics{i};
    [acc, loss] = CBIG_compute_prediction_acc_and_loss(pred,y_test, curr_metric,y_train);
    optimal_stats.(curr_metric) = acc;
end

acc = optimal_stats.corr;

%% testing
fprintf('>>> Test accuracy (correlation): %f\n\n', acc);

end