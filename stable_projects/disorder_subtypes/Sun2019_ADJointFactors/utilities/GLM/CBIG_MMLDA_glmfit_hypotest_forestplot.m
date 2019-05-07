function p_pairwise = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, out_name, compare_factor_names, CI_scale, p_thresh)
% CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, out_name, compare_factor_names, CI_scale, p_thresh)
%
% This function will do glmfit, hypothesis tests and forest plot.
%
% Input:
%   - k             : number of factors
%   - X             : regressors without constant columns 1s
%   - y             : target variable 
%   - dependency    : if 'dependent' which means the factors sum to one. 
%                     the X should include only include F2, F3 not F1.
%                     if 'independent' the factors are independent to each other 
%                     the X should include F1, F2, F3
%   - out_dir       : out_name directory
%   - compare_factor_names : cell array. Comparison between factors, such as {'F2-F1'}
%   - CI_scale      : (optional) confidence interval scale. Default 1.96 (95% CI)
%   - p_thresh      : (optional) FDR p value threshold. Default 0.05. 
%
% Example:
%   CBIG_MMLDA_glmfit_hypotest_forestplot(3, X, y, 'dependent', './', 'age', {'F2-F1'})
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 9
	p_thresh = 0.05;
end
if nargin < 8
	CI_scale = 1.96;
end

mkdir(out_dir)

% exclude NaN and Inf in y
y_nan = y(~(isnan(y)|isinf(y)));
X_nan = X(~(isnan(y)|isinf(y)),:);

N = length(y_nan);

% out_name to console and txt file
disp(['Length of y excluding NaN and Inf:' num2str(length(y_nan))])

% fit GLM model
[b, ~, stats] = glmfit(X_nan, y_nan);

% hypo test
if strcmp(dependency, 'independent')
	result = CBIG_MMLDA_hypotest_independent_factors(k, N, b, stats, [out_dir '/' out_name '.txt'], CI_scale)
elseif strcmp(dependency, 'dependent')
	result = CBIG_MMLDA_hypotest_dependent_factors(k, N, b, stats, [out_dir '/' out_name '.txt'], CI_scale)
else
	error('No such option')
end
p_pairwise = result(:, end);

% plot the forest figure 
CBIG_MMLDA_forest_plot(result, 'auto-center', [out_dir '/' out_name]);

% add text to forest figure 
CBIG_MMLDA_add_text_to_forest_plot([out_dir '/' out_name '.png'], [out_dir '/' out_name '.txt'], [out_dir '/'], p_thresh, k, compare_factor_names)