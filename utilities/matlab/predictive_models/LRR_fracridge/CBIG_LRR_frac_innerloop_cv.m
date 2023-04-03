function [loss, y_pred, acc_corr_fold] = CBIG_LRR_frac_innerloop_cv(feature_mat, y, lambda, tot_folds, metric)

% [loss_mean] = CBIG_LRR_frac_innerloop_cv(feature_mat, y, lambda, tot_folds, metric)
% 
% This function performs the inner-loop cross-validation of linear ridge
% regression using the fitridge function provided by:
% https://github.com/nrdg/fracridge. This function is an adaptation of 
% CBIG_LRR_innerloop_cv.
% 
% Inputs:
%   - feature_mat:
%     A 2D matrix of the selected features (dim: #subjects x #features).
%     This function automatically prepends the bias term (a vector of ones)
%     onto "feature_mat".
% 
%   - y:
%     A matrix of the target variables to be predicted (dimension:
%     #subjects x K). K is the number of target variables to be predicted.
% 
%   - lambda:
%     A vector, the fraction hyperparameter. Value is between 0 and 1.
% 
%   - tot_folds:
%     A scalar, how many inner-loop CV folds the user wants to conduct.
%
%   - metric
%     A string, what is the accuracy metric that we want to optimise our
%     hyperparameter with. Can choose from 'corr','COD','predictive_COD',
%     'MAE_norm','MAE','MSE','MSE_norm'. Default is 'predictive_COD'.
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
% Written by Yanrui, Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
        'non_default_packages', 'fracridge'));

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
    
    y_train = y(curr_train_idx,:);
    y_test = y(curr_test_idx,:);
    
    % training using fractional ridge regression with specific lambda value
    
    beta = fracridge(feat_train,lambda,y_train); %beta is a #feature x #lambda x #behaviors
    beta = permute(beta, [1 3 2]); %beta is a #feature x #behaviors x #lambda 
    % prediction
    for l = 1:length(lambda)
        pred = feat_test * beta(:,:,l); %pred: #sub x #behaviors
        y_pred_tmp(:,:,l) = pred;
        

        for i = 1:size(pred, 2)
            [acc_corr_fold(curr_fold,i,l), loss(curr_fold,i,l)] = CBIG_compute_prediction_acc_and_loss(pred(:,i),...
            y_test(:,i), metric,y_train(:,i));
        end
    end
    y_pred{curr_fold} = y_pred_tmp;
    clear y_pred_tmp
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
'non_default_packages', 'fracridge'));
end

