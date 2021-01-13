function [loss_mean, y_pred] = CBIG_LRR_fitrlinear_innerloop_cv( feature_mat, y, lambda, tot_folds, metric )

% [loss_mean] = CBIG_LRR_fitrlinear_innerloop_cv( feature, y, lambda, tot_folds, metric )
% 
% This function performs the inner-loop cross-validation of linear ridge
% regression using the fitrlinear function provided by MATLAB. This function is an
% adaptation of CBIG_LRR_innerloop_cv.
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
%   - metric
%     A string, what is the accuracy metric that we want to optimise our
%     hyperparameter with. Can choose from
%     'corr','COD','predictive_COD','MAE_norm','MAE','MSE','MSE_norm'
% 
% Outputs:
%   - loss_mean:
%     A scalar. It is the mean negative loss from the given metric averaged 
%     across inner-loop folds.
%
%  - y_pred:
%    A cell of predicted y values from each fold. There are #tot_folds cells with
%    #subjects / #tot_folds predictions in each cell.
% 
% Written by Yanrui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(size(y,2)~=1)
    error('''y''should be a column vector.')
end


%% split data
rng('default');  rng(1);
num_sub = size(feature_mat, 1);
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
    
    % training using fitrlinear with specific lambda value , least squares
    % and ridge regression
    
    Mdl = fitrlinear(feat_train',y_train,'ObservationsIn','columns', 'Lambda',lambda, 'Learner',...
        'leastsquares', 'Regularization','ridge');
   
    % prediction
    pred = predict(Mdl,feat_test','ObservationsIn','columns');
    y_pred{curr_fold} = pred;
    
    [acc_corr_fold(curr_fold,1), loss(curr_fold,1)] = CBIG_compute_prediction_acc_and_loss(pred,y_test,...
        metric,y_train);
    
end

loss_mean = mean(loss);
loss_mean = -loss_mean;
fprintf('Mean testing loss (%s) across %d folds: %f\n', metric, tot_folds, loss_mean);

end

