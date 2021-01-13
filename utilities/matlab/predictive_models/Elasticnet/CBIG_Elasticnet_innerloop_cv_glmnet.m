function [mean_acc, neg_loss, lambdas, y_pred_all] = ...
CBIG_Elasticnet_innerloop_cv_glmnet( feature_mat, y, alpha, lambda_set, tot_folds, metric )

% This function performs the inner-loop cross-validation of generalised linear 
% model using the glmnet function provided in `http://web.stanford.edu/~hastie/glmnet_matlab/`. 
% This function is an adaptation of CBIG_LRR_innerloop_cv.
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
%   - alpha: 
%     A scalar, the ratio of L1 to L2 regularisation
%
%   - lambda_set:
%     A vector of regularization hyperparameters for the grid search.
%     If the set passed in is empty, glmnet will initialise 200 
%     lambda values for the gridsearch.
%
%   - tot_folds:
%     A scalar, how many inner-loop CV folds the user wants to conduct.
%
%   - metric:
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
% 
% Outputs:
%   - mean_acc:
%     A scalar. It is the mean accuracy averaged across
%     inner-loop folds in terms of the chosen input metric.
%
%   - neg_loss:
%     A scalar. It is the mean negative loss averaged across
%     inner-loop folds in terms of the chosen input metric.
%     Higher values indicate greater accuracy.
%
%   - lambdas:
%     A vector. It is the final set of lambdas used for the
%     inner loop. If lambda_set was passed in, it will be 
%     the same as lambda_set.
%
%   - y_pred_all:
%     A cell containing the final predicted target variable 
%     for each fold of the inner loop. 
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(size(y,2)~=1)
    error('''y''should be a column vector.')
end


%% split data
rng('default');  rng(1);
num_sub = size(feature_mat, 1);
p = size(feature_mat, 2);
cv_idx = cvpartition(num_sub, 'kfold', tot_folds);

foldid = zeros(size(y));

for i = 1:tot_folds
    foldid(cv_idx.test(i),:) = i;
end

%% For each fold, training and testing
 
% training
for curr_fold = 1:tot_folds
    
    curr_train_idx = cv_idx.training(curr_fold);
    curr_test_idx = cv_idx.test(curr_fold);
    
    feat_train = feature_mat(curr_train_idx, :);
    feat_test = feature_mat(curr_test_idx, :);
    
    y_train = y(curr_train_idx);
    y_test = y(curr_test_idx);
    
    %% update elasticnet parameters
    if ~isempty(lambda_set) % use preset range of lambda
      opts.lambda = lambda_set;
      opts.alpha = alpha;    
      fit=glmnet(feat_train,y_train, [], opts);    
      y_predict = glmnetPredict(fit, feat_test); 
      lambdas = lambda_set;
    else %% let glmnet choose the lambda
      opts.alpha = alpha;
      if curr_fold == 1    
         opts.nlambda = 200; % set lambdas to first 200 values or less
         fit=glmnet(feat_train,y_train, [], opts);  
         lambdas = fit.lambda;
      else
         opts.lambda = lambdas;
         fit=glmnet(feat_train,y_train, [], opts);
      end
      y_predict = glmnetPredict(fit, feat_test); 
    end
         
     % compute accuracy
     [acc_fold(curr_fold, :), loss_fold(curr_fold, :)] = CBIG_compute_prediction_acc_and_loss...
          (y_predict, y_test, metric, y_train);
     y_pred_all{curr_fold} = y_predict;
    
end
    
% return mean cv accuracy for all lambdas
mean_acc = mean(acc_fold); 
neg_loss = -mean(loss_fold);

end

