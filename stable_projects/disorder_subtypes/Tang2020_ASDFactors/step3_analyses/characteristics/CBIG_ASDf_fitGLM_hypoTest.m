function num_pVals = CBIG_ASDf_fitGLM_hypoTest(k, X, y, numReg, output_name)
% num_pVals = CBIG_ASDf_fitGLM_hypoTest(k, X, y, numReg, output_name)
%
% This function fits GLM and performs hypothesis test for the GLM, and
% plots the results. If output_name is specified, the plot will be saved.
%
% Input:
%     - k:
%           Integer. Number of factors
%     - X:
%           NxA matrix, where N is the number of subjects and A is the number 
%           of explanatory variables + nuisance variables
%     - y:
%           Nx1 matrix, where N is the number of subjects
%     - numReg:
%           Number of regressors
%     - output_name (optional):
%           File name (including path) that will be used to save the plot
% Output:
%     - num_pVals:
%           (V+1)x1 vector, where V is the number of pairwise comparisons.
%           E.g., if k=3, V is 3; if k=4, V is 6. 1st row stores the number
%           of subjects in this logistic regression analysis, 2nd to last
%           rows store the pairwise p-values.
%
% Example:
%     num_pVals = CBIG_ASDf_fitGLM_hypoTest(3, X, y, 13, 'age_glm')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if size(y,2) > 1
    error('Input argument "y" should be a column vector.');
end

if size(X,1) ~= size(y,1)
    error('Input argument "X" should have dimension N x A, where N is the number of subjects.');
end

no_comparisons = nchoosek(k,2);
num_pVals = zeros(no_comparisons+1,1);

%% exclude NaN and Inf in y
y_nan = y(~(isnan(y)|isinf(y)));
X_nan = X(~(isnan(y)|isinf(y)),:);

% Display in console
disp(['Length of y excluding NaN and Inf:' num2str(length(y_nan))])
num_pVals(1) = length(y_nan);

%% Fit GLM model and perform hypothesis test
[b, ~, stat] = glmfit(X_nan, y_nan);
% stat.dfe
size(stat.covb)
[result, num_pVals(2:end,1)] = CBIG_ASDf_hypoTest(b, stat, numReg, k);

%% Plot the results
if nargin > 4 && ~isempty(output_name)
    CBIG_ASDf_plotCmp(result, 'auto-center', 0, output_name);
else
    CBIG_ASDf_plotCmp(result, 'auto-center', 0);
end
