function [acc_corr, mse, y_predict] = CBIG_LRR_train_test( feat_train, feat_test, y_train, y_test, lambda )

% [acc_corr, mse, y_predict] = CBIG_LRR_train_test( feat_train, feat_test, y_train, y_test, lambda )
% 
% Based on the specified regularization hyperparameter (lambda), this
% function trains a linear ridge regression model to predict the target
% value (y_train) from the features (feat_train) in the training set. After
% that, it applies the trained model on the test set (y_test and feat_test).
% 
% Inputs:
%   - feat_train:
%     A 2D matrix (#training subjects x #features). It is the selected
%     features of training subjects. A column vector of all-ones will be
%     automatically prepended to feat_train to capture the interception
%     term.
% 
%   - feat_test:
%     A 2D matrix (#test subject x x#features). It is the selected features
%     of test subjects. A column vector of all-ones will be automatically
%     prepended to feat_test to capture the interception term.
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
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(size(y_train,2) > 1)
    error('''y_train'' should be a column vector')
end
if(size(y_test,2) > 1)
    error('''y_test'' should be a column vector')
end

feat_train = cat(2, ones(size(feat_train,1),1), feat_train);
feat_test = cat(2, ones(size(feat_test,1),1), feat_test);

p = size(feat_train, 2);

%% training
beta = (feat_train' * feat_train + lambda*eye(p))' \ feat_train' * y_train;
y_predict_train = feat_train * beta;
acc_corr_train = CBIG_corr(y_train, y_predict_train);
mse_train = mean((y_train - y_predict_train).^2);
fprintf('>>> Training accuracy (correlation): %f, mse: %f\n', acc_corr_train, mse_train);

%% testing
y_predict = feat_test * beta;
acc_corr = CBIG_corr(y_test, y_predict);
mse = mean((y_test - y_predict).^2);
fprintf('>>> Test accuracy (correlation): %f, mse: %f\n\n', acc_corr, mse);


end

