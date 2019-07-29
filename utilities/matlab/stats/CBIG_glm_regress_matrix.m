function [resid_mtx, coef_mtx, std_mtx, retrend_mtx] = CBIG_glm_regress_matrix(input_mtx, regressor, polynomial_fit, censor)

% [resid_mtx, coef_mtx, std_mtx, retrend_mtx] = CBIG_glm_regress_matrix(input_mtx, regressor, polynomial_fit, censor)
% 
% This function does GLM regression for each column of a MxN input matrix 
% (input_mtx) independently with a MxK regressor matrix (regressor). If
% polynomial_fit is set to -1, we add nothing to the regressor matrix. If 
% polynomial_fit is set to 0, we prepend a Mx1 vector [1,1,1...1]' to the regressor
% matrix and end up with a Mx(K+1) matrix. If polynomial_fit is set to 1, 
% a Mx1 vector [1,1,1...1]' will be added into the first column and a Mx1 
% matrix linspace(-1, 1, M)' will be added into the second column of the regressor
% matrix. 
% If censor is given, this function will first exclude the censored frames 
% from input matrix and regressor matrix: input_mtx(censor, :), 
% regressor(censor, :), then calculate the coefficients (coef_mtx) and 
% apply it to original input matrix (input_mtx) and regressor matrix 
% (regressor) to get residual matrix (resid_mtx). 
% 
% Input:
%     - input_mtx:
%       M x N matrix, where M is dimension of a sample, N is number of
%       sample.
%
%     - regressor:
%       M x K matrix, where M is dimension of a sample, K is number of
%       regressors.
%
%     - polynomial_fit:
%       -1/0/1/ (default is 0). If polynomial_fit is set to -1, we add
%       nothing in regressor matrix. If polynomial_fit is set to 0, we prepend a 
%       Mx1 vector [1,1,1...1]' to the regressor matrix and end up with a 
%       Mx(K+1) matrix. If polynomial_fit is set to 1, a Mx1 vector 
%       [1,1,1...1]' will be added into the first column and a Mx1 matrix 
%       linspace(-1, 1, M)' will be added into the second column of the 
%       regressor matrix.
%
%     - censor:
%       a Mx1 vector with 0 and 1 (default is []), 1 means kept frames, 0 
%       means removed frames. If the censor vector is empty, this step will
%       not remove any frames. If the censor vector is given, this function 
%       will remove the censored frames from the input and regressor matrix.
%
% Output:
%     - resid_mtx:
%       a M x N matrix, where M is dimension of a sample, N is number of 
%       samples. Resulted residual matrix after GLM regression.
% 
%     - coef_mtx:
%       a (K+1) x N matrix, if linear_detrend is set to 0, the first row
%       of coefficient matrix corresponds to regressor [1,1,1...1]'. 
%       a (K+2) x N matrix, if linear_detrend is set to 1, the second row 
%       of coefficient matrix corresponds to regressor linspace(-1, 1, M)'.
%       K is number of input regressors and N is number of samples. 
%       Resulted Coefficient matrix after GLM regression.
% 
%     - std_mtx:
%       a 1 x N matrix, where N is number of samples. Standard deviation
%       of residual matrix.
% 
%     - retrend_mtx:
%       a M x N matrix, where M is dimension of a sample, N is number of
%       samples. If input regressor is empty and polynomial is set to 3, 
%       this function actually demean and remove the linear trend of input
%       matrix. retrend_mtx = Mx2 regressor matrix X 2xN coefficent matrix.
%       If someone wants to add back the linear trend of original 
%       data, the retrend_mtx can be used to add to the residual matrix. 
%       If input regressor is not empty, retrend is not allowed and
%       retrend_mtx is empty.       
%  
% Example:
% [resid_mtx, coef_mtx, std_mtx, retrend_mtx] = CBIG_glm_regress_matrix(input_data_mtx, regressor_mtx, 1)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% default polynomial_fit option is 1,  censor option is []

if (nargin <3)
    polynomial_fit = 0;
end
if (nargin < 4)
    censor = [];
end

if ~isempty(censor)
    if size(censor,2) ~= 1
        error('Input argument ''censor'' should be a column vector');
    end
end

% check if there are no regressor
if ((isempty(regressor)) && (polynomial_fit == -1))
    error('ERROR: No regressor, quit!')
end

%% construct GLM regressors 
% If polynomial_fit is set to 0, we will prepend a Mx1 vector 
% [1,1,1...1]' to the regressor matrix. X is resulted Mx(K+1) matrix.
% If polynomial_fit is set to 1, a Mx1 vector [1,1,1...1]' will be added 
% into the first column and a Mx1 matrix linspace(-1, 1, M)' will be added 
% into the second column of the regressor matrix.
Y = input_mtx;
X = regressor;
if (polynomial_fit == 1)
    X = [linspace(-1, 1, size(Y, 1))' X]; 
end
if ((polynomial_fit == 0) || (polynomial_fit == 1))
    X = [ones(size(Y, 1), 1) X];
end

%% do least square GLM regression for input matrix Y, caculate coefficient matrix b
if (~isempty(censor))
    % if censor vector is not empty
    censor = logical(censor);
    Y_censor = Y(censor, :);
    X_censor = X(censor, :);
    b = double(X_censor'*X_censor)\double(X_censor'*Y_censor);
else
    % if censor vector is empty
    b = double(X'*X)\double(X'*Y);
    %b = (X'*X)\(X'*Y);
end

%% Output coefficient matrix, residual matrix and retrend_matrix
resid_mtx = Y - X*b;
coef_mtx = b;
std_mtx = std(resid_mtx, 0, 1);
if (isempty(regressor))
    % if regressor is empty
    retrend_mtx = X*b;
else
    % if regressor is not empty
    retrend_mtx = [];
end
