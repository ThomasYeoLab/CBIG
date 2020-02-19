function [logReg_wMeanStd, num_pVal] = CBIG_ASDf_logReg_compareSex(k, sub_id, ...
sub_info_file, factorLoading_file, factor_order, output_name)
% [logReg_wMeanStd, num_pVal] = CBIG_ASDf_logReg_compareSex(k, sub_id, 
% sub_info_file, factorLoading_file, factor_order, output_name)
% 
% This function performs logistic regression to compare ASD participants' sex
% across latent factors, and plot the errorbar plot. If output_name is
% specified, the plot will be saved.
% NOTE: This function uses some functions in 
% $CBIG_CODE_DIR/stable_projects/disorder_subtypes/Zhang2016_ADFactors.
%
% Input:
% 	  - k:
%           Integer. Number of factors
%     - sub_id:
%           Nx1 cell array, where N is the number of subjects
%     - sub_info_file:
%           .csv file containing all subjects' demographics info
%     - factorLoading_file:
%           .txt file containing all ASD sujects' factor compositions
%     - factor_order:
%           The order of the factors. e.g., [1 3 2] means re-ordering the
%           factors, 1st factor as factor 1, 3rd factor as factor 2, and 2nd
%           factor as factor 3.
%     - output_name (optional):
%           File name (including full path) to save the plot
% Output:
%     - logReg_meanStd:
%           kx1 cell array, where k is the number of factors. The i-th entry
%           is the weighted mean and weighted standard
%           deviation of the i-th factor.
%     - num_pVal:
%           (V+1)x1 vector, where V is the number of pairwise comparisons.
%           E.g., if k=3, V is 3; if k=4, V is 6. 1st row stores the number
%           of subjects in this logistic regression analysis, 2nd to last
%           rows store the pairwise p-values.
% 
% Example:
%       [logReg_wMeanStd, num_pVal] = CBIG_ASDf_logReg_compareSex(3, sub_id, 
%       '~/example_input/subInfo.csv', '~/example_output/visualizeFactors/factorComp.txt', [3 2 1], 'k3_sex')
% 
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variable
if size(factor_order,1) > 1
    error('Input argument "factor_order" should be a row vector.');
end

logReg_wMeanStd = cell(k,1);
no_comparisons = nchoosek(k,2);
num_pVal = zeros(no_comparisons+1,1);

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ZHANG_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Zhang2016_ADFactors');
addpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','functions'));
addpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','characteristics'));
addpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','MCI2ADProgression'));

%% Get sex & factor loadings
[~, ~, ~, id_sex, ~, ~, ~, id_factorLoading] = CBIG_ASDf_getSubData(sub_info_file, sub_id, ... 
[], [], factorLoading_file);

sex = cell2mat(id_sex(:,2));
factorComp = cell2mat(id_factorLoading(:,2:end));
factorComp = factorComp(:,factor_order);

%% Get regressors
regressors = CBIG_ASDf_genRegressors(sub_info_file, sub_id);
site_reg = regressors(:,4:end); % Sites as regressors

%% Convert male to 0, female to 1
y = sex - 1; % Male: 0; Female: 1

%% Run logistic regression
[wMean, wStdDev] = CBIG_compute_weightedMeanStdDev(factorComp, y); % use Xiuming's function

num_pVal(1) = length(y);

%%% K=2
if k == 2
    % Hypothesis testing: exact LR test
    X = [factorComp(:, 1) site_reg];
    n = ones(size(y, 1), 1);
    [b, ~, stats] = glmfit(X, [y n], 'binomial', 'link', 'logit');
    disp('Exact LR test');
    % Unrestricted model
    b_u = b;
    yFit = glmval(b_u, X, 'logit', 'size', n);
    ll_u = sum(log(binopdf(y, n, yFit./n)));
    
    % Forest plot data
    A = zeros(1,2);
    % F1 - F2
    H = [0 1];
    idx1 = find(H==1);
    A(1,:) = [-stats.beta(idx1), sqrt(stats.covb(idx1,idx1))];
    
    % Overall
    % Restricted model
    X_r = site_reg;
    dof = abs(size(X,2) - size(X_r,2));
    disp('Overall');
    p = CBIG_lr_test(X_r, y, n, ll_u, dof);
    fprintf('p = %e\n', p);
    num_pVal(2) = p;
    fprintf('Factor1 mean(std): %.3f (%.3f)\n', wMean(1), wStdDev(1));
    fprintf('Factor2 mean(std): %.3f (%.3f)\n', wMean(2), wStdDev(2));
    logReg_wMeanStd{1} = sprintf('Factor1 mean(std): %.3f (%.3f)\n', wMean(1), wStdDev(1));
    logReg_wMeanStd{2} = sprintf('Factor2 mean(std): %.3f (%.3f)\n', wMean(2), wStdDev(2));
    
%%% K=3
elseif k == 3
    % Hypothesis testing: exact LR test
    X = [factorComp(:, [1 2]) site_reg];
    n = ones(size(y, 1), 1);
    [b, ~, stats] = glmfit(X, [y n], 'binomial', 'link', 'logit');
    disp('Exact LR test');
    
    % Unrestricted model
    b_u = b;
    yFit = glmval(b_u, X, 'logit', 'size', n);
    ll_u = sum(log(binopdf(y, n, yFit./n)));
    
    % Forest plot data
    A = zeros(3,2);
    % F1 - F3
    H = [0 1 0];
    idx1 = find(H==1);
    A(1,:) = [-stats.beta(idx1), sqrt(stats.covb(idx1, idx1))];  
    % F2 - F3
    H = [0 0 1];
    idx1 = find(H==1);
    A(2,:) = [-stats.beta(idx1), sqrt(stats.covb(idx1, idx1))]; 
    % F2 - F1
    H = [0 -1 1];
    idx1 = find(H==1);
    idx2 = find(H==-1);
    [mu, sigma2] = CBIG_compute_mu_sigma2(stats, idx1, idx2);
    A(3,:) = [-mu, sqrt(sigma2)];
    
    % Overall
    % Restricted model
    X_r = site_reg;
    dof = size(X,2) - size(X_r,2);
    disp('Overall');
    p = CBIG_lr_test(X_r, y, n, ll_u, dof);
    fprintf('p = %e\n', p);
    % Restricted model
    X_r = [X(:,2) site_reg];
    dof = size(X,2) - size(X_r,2);
    disp('Factor1 == Factor3?');
    p = CBIG_lr_test(X_r, y, n, ll_u, dof);
    fprintf('p = %e\n', p);
    num_pVal(2) = p;
    % Restricted model
    X_r = [X(:,1) site_reg];
    dof = size(X,2) - size(X_r,2);
    disp('Factor2 == Factor3?');
    p = CBIG_lr_test(X_r, y, n, ll_u, dof);
    fprintf('p = %e\n', p);
    num_pVal(3) = p;
    % Restricted model
    X_r = [X(:,1) + X(:,2) site_reg];
    dof = size(X,2) - size(X_r,2);
    disp('Factor1 == Factor2?');
    p = CBIG_lr_test(X_r, y, n, ll_u, dof);
    fprintf('p = %e\n', p);
    num_pVal(4) = p;

    fprintf('Factor1 mean(std): %.3f (%.3f)\n', wMean(1), wStdDev(1));
    fprintf('Factor2 mean(std): %.3f (%.3f)\n', wMean(2), wStdDev(2));
    fprintf('Factor3 mean(std): %.3f (%.3f)\n', wMean(3), wStdDev(3));
    logReg_wMeanStd{1} = sprintf('Factor1 mean(std): %.3f (%.3f)\n', wMean(1), wStdDev(1));
    logReg_wMeanStd{2} = sprintf('Factor2 mean(std): %.3f (%.3f)\n', wMean(2), wStdDev(2));
    logReg_wMeanStd{3} = sprintf('Factor3 mean(std): %.3f (%.3f)\n', wMean(3), wStdDev(3));
else
    fprintf('NOT CONFIGURED! Choose from K = 2 or K = 3.');  
end

if nargin > 5 && ~isempty(output_name)
    CBIG_ASDf_plotCmp(A, 'auto-center', 1, output_name);
else
    CBIG_ASDf_plotCmp(A, 'auto-center', 1);
end

%% Remove paths
rmpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','functions'));
rmpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','characteristics'));
rmpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','MCI2ADProgression'));
