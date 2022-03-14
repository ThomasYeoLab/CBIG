function y_hat = CBIG_regress_X_from_y_test(y, X, beta, X_train_mean)

% y_hat = CBIG_regress_X_from_y_test(y, X, beta, X_train_mean)
% 
% This function applies the regression coefficients "beta" generated from
% the training set (see CBIG_regress_X_from_y_train.m) to regress the
% nuisance variables X from the target variables y.
% 
% Inputs:
%   - y
%     A #TestSubjects x 1 vector of target variables in the test set.
% 
%   - X
%     A #TestSubjects x #covariates of regressors in the test set. The
%     regressors will be automatically demeaned. A column vector of ones
%     will be automatically prepended for the bias term.
%     X can be empty. In that case, this function only regresses a vector
%     of ones from y (i.e. demean).
% 
%   - beta
%     A #covariates x 1 vector of regression coefficients used to regress
%     the training set target variables. It is derived from the output of
%     "CBIG_regress_X_from_y_train".
%
%   - X_train_mean
%     A #covariates x 1 vector of the mean value of each covariate calculated
%     from the training set. 
% 
% Outputs:
%   - y_hat
%     A #TestSubjects x 1 vector of residual target variable values after
%     regression in the test set.
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Authors: Ru(by) Kong,  Jingwei Li

% check y in case there is nan
nan_index_y = isnan(sum(y,2));
% check X in case there is nan, exclude this subject for all covariates
nan_index_X = isnan(sum(X,2));

if(~isempty(X))
    nan_index = (nan_index_y | nan_index_X);
else
    nan_index = nan_index_y;
end

% keep original y so that y has the same length for all measures
y_hat = y;
y_hat(nan_index,:) = NaN;

y(nan_index,:) = [];
if(~isempty(X))
    X(nan_index,:) = [];
end

% demean nusiance regressors
X = bsxfun(@minus, X, X_train_mean);

% Add bias term for X

X = [ones(size(y,1),1) X];
% X = ones(size(X,1),1);

% y = X*beta + e
% 
y_hat(~nan_index,:) = y - X*beta;

end

