function [N, wMean, wStdDev, p] = CBIG_MMLDA_characteristics_glm(prob, y, y_type)
% [wMean, wStdDev, p] = CBIG_MMLDA_characteristics_glm(prob, y, y_type) 
%
% Analysis about factors effect on characteristics by using glmfit and hypothesis test.
%
% Input:
%   - prob      : 
%           N x K matrix, factor compositions. N is number of subjects, K is number of factors. 
%   - y         : 
%           N x 1 vector, target characteristics. 
%   - y_type (optional) : 
%           if 'continuous' (default), will use general linear model and linear hypothesis test.
%           if 'binary', will use logistic regression model and likelihood ratio test. 
%
% Output:
%   - N         :
%           number of subjects after removing missing characteristics
%   - wMean     :
%           1 x K vector, weighted mean of characteristics for each factor
%   - wStdDev   :
%           1 x K vector, weighted standard deviation of characteristics for each factor
%   - p         :
%           scalar, p value for omnibus test
%
% Example:
%   [wMean, wStdDev, p] = CBIG_MMLDA_characteristics_glm(prob, y, 'continuous')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if nargin < 3
    y_type = 'continuous';
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step3_analyses_internalUse/characteristics'])
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step3_analyses_internalUse/functions'])

% exclude NaN and Inf in y
y_nan = y(~(isnan(y)|isinf(y)));
prob_nan = prob(~(isnan(y)|isinf(y)),:);

N = length(y_nan);

% compute weighted mean and std
[wMean, wStdDev] = CBIG_compute_weightedMeanStdDev(prob_nan, y_nan);

% fit GLM model
k = size(prob_nan, 2);
X = prob_nan(:, 2:end);

if strcmp(y_type, 'continuous')
    [b, ~, stats] = glmfit(X, y_nan);

    % Overall (F2 - F1 = 0; F3 - F1 = 0)
    H = [0 1 0; 0 0 1];
    c = [0; 0]; % H*b = c
    p = linhyptest(b, stats.covb, c, H, stats.dfe);
    % Pairwise comparisons
    % F2 - F1
    H = [0 1 0];
    c = 0;
    if H*b > 0
        fprintf('F2 > F1\n');
    else
        fprintf('F2 < F1\n');
    end
    p1 = linhyptest(b, stats.covb, c, H, stats.dfe);
    fprintf('p = %e\n', p1);
    
    % F3 - F1
    H = [0 0 1];
    c = 0;
    if H*b > 0
        fprintf('F3 > F1\n');
    else
        fprintf('F3 < F1\n');
    end
    p2 = linhyptest(b, stats.covb, c, H, stats.dfe);
    fprintf('p = %e\n', p2);
    
    % F3 - F2
    H = [0 -1 1];
    c = 0;
    if H*b > 0
        fprintf('F3 > F2\n');
    else
        fprintf('F3 < F2\n');
    end
    p3 = linhyptest(b, stats.covb, c, H, stats.dfe);
    fprintf('p = %e\n', p3);
elseif strcmp(y_type, 'binary')
    n = ones(size(y_nan, 1), 1);
    [b, ~, stats] = glmfit(X, [y_nan n], 'binomial', 'link', 'logit'); % log(�/(1-�)) = Xb
    % disp('Exact LR test');
    % Unrestricted model
    b_u = b;
    yFit = glmval(b_u, X, 'logit', 'size', n);
    ll_u = sum(log(binopdf(y_nan, n, yFit./n)));
    % Overall
    % Restricted model
    X_r = [];
    dof = 2;
    % disp('Overall');
    p = CBIG_lr_test(X_r, y_nan, n, ll_u, dof);
end

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step3_analyses_internalUse/characteristics'])
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step3_analyses_internalUse/functions'])
