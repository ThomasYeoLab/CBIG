function [y_hat, beta] = CBIG_regress_X_from_y_train(y, X)

% [y_hat, beta] = CBIG_regress_X_from_y_train(y, X)
% 
% This function regress covariates X (e.g. head motion) from y (e.g.
% behavioral scores) in the training set (for any preditive model).
% 
% Inputs:
%   - y 
%     A #TrainingSubjects x 1 vector.  Each row in the matrix contains the
%     target variable value of each subject.
% 
%   - X
%     A #TrainingSubjects x #covariates matrix of regressors in the training set.
%     Each row in the matrix contains the regressor values for each subject.
%     The regressors will be automatically demeaned. A column vector of
%     ones will be prepended automatically for the bias term.
% 
% Outputs:
%   - y_hat
%     A #TrainingSubjects x 1 vector of residual target variable values
%     after regression.
% 
%   - beta
%     A #covariates x 1 vector of regression coefficients used to regress
%     the training set target variables.
%
% Written by Ru(by) Kong, Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% check y in case there is nan
nan_index_y = isnan(y);
% check X in case there is nan, exclude this subject for all covariates
nan_index_X = isnan(sum(X,2));

nan_index = (nan_index_y | nan_index_X);

% keep original y so that y has the same length for all measures
y_hat = y;
y_hat(nan_index,:) = NaN;

y(nan_index) = [];
X(nan_index,:) = [];

% demean nusiance regressors
X = bsxfun(@minus, X, mean(X));

% Add bias term for X

X = [ones(size(X,1),1) X];
% X = ones(size(X,1),1);

% y = X*beta + e
% 
% estimate beta
%  beta = ( X' * X )^-1 * X' * y
%  Px1  =  PxN  NxP      PxN  Nx1
beta = (X'*X)\(X'*y);
y_hat(~isnan(y_hat)) = y - X*beta;

end