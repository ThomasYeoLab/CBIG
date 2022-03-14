function [acc_corr, y_predict, y_pred_train, optimal_stats] = ...
    CBIG_Elasticnet_train_test_glmnet( feat_train, feat_test, y_train, y_test, ...
    alpha, lambda)
 
% Based on the specified regularization ratio (alpha) and regularisation hyperparameter (lambda), 
% this function trains a generalized linear model to predict the target value (y_train)
% from the features (feat_train) in the training set. After that, it applies the trained model
% on the test set (y_test and feat_test). This function is an adaptation from CBIG_LRR_train_test
% in CBIG_repo. It uses the MATLAB function `glmnet` from  `http://web.stanford.edu/~hastie/glmnet_matlab/`
% for fitting the generalized linear model. The glmnet function automatically fits a bias term
% from input features. Hence, there is no need for the user to concatenate
% a bias term manually.
%
% Predictions are made for the following statistics:
%
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
%   - alpha: 
%     A scalar, the ratio of L1 to L2 regularisation
%
%   - lambda:
%     A scalar, the regularization hyperparameter.
%    
% 
% Outputs:
%   - acc_corr:
%     A scalar. It is the correlation accuracy in the test set.
%   - y_predict:
%     A vector. It is the predicted values of the test set using the model
%     fitted from the training set.
%   - y_pred_train:
%     A vector. It is the predicted values of the training set using the model
%     fitted from the training set. It can be used for model
%     interpretation.
%   - optimal_stats:
%     A struct with 7 fields (one for each test statistic, see above). 
%     Each field contans a scalar of the test statistic accuracy.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(size(y_train,2) > 1)
    error('''y_train'' should be a column vector')
end
if(size(y_test,2) > 1)
    error('''y_test'' should be a column vector')
end

p = size(feat_train, 2);
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};


%% training
rng('default')
rng(1)
opts.lambda = lambda; % should take best lambda value from innerloop cv here
opts.alpha = alpha;
fit=glmnet(feat_train,y_train, [], opts);

%% testing
y_predict = glmnetPredict(fit, feat_test);
y_pred_train = glmnetPredict(fit, feat_train);

% Print warning standard deviation of prediction is close to 0
if std(y_predict) < 2e-10
    fprintf('>>> WARNING: Predictions are similar across subjects, consider lowering lambda \n')
end

for i = 1:length(metrics)
    curr_metric = metrics{i};
    [acc, loss] = CBIG_compute_prediction_acc_and_loss(y_predict ,y_test, curr_metric,y_train);
    optimal_stats.(curr_metric) = acc;
end

acc_corr = optimal_stats.corr;

fprintf('>>> Test accuracy (corr): %f\n\n', acc_corr);


end

