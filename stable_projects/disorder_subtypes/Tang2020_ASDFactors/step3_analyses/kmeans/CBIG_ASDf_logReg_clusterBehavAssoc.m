function [model, p_lr] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, ...
sub_info_file, behav_score_file, idx_score)
% [model, p_lr] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, 
% sub_info_file, behav_score_file, idx_score)
%
% Performs logistic regression to compare ASD participants' behavioral
% scores across k-means clusters, and plot the errorbar plot. If output_name is
% specified, the plot will be saved.
%
% Input:
% 	  - cluster_idx:
%           Integer, index of cluster, e.g. 1, 2 or 3
%     - id_scoreGrp:
%           IDs of subjects of interests. Nx1 cell array, where N is the
%           number of subjects
%     - sub_info_file:
%           .csv file containing all subjects' demographics info
%     - behav_score_file:
%           .csv file containing all subjects' behavioral scores
%     - idx_score:
%           Column indices of behavioral scores of interests in
%           behav_score_file
%
% Output:
%     - model:
%           Fitted GLM model from fitglm function
%     - p_lr:
%           p-values of model statistical significance for all logistic regression models
%
% Example:
%       [model, p_lr] = CBIG_ASDf_logReg_clusterBehavAssoc(1, id_scoreGrp, '~/input/subInfo.csv', 
%                       '~/input/behav_scores.csv', [1 2 3])
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variable
if size(idx_score,1) > 1 && size(idx_score,2) > 1
    error('Input argument "idx_score" should be a row or column vector.');
end

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ZHANG_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Zhang2016_ADFactors');
addpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','functions'));
addpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','characteristics'));
addpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','MCI2ADProgression'));

%% Get participants' IDs in the current cluster having the set of behavior scores
[~, id_dx] = CBIG_ASDf_getSubData(sub_info_file);
id_all = id_dx(:,1);
dx = cell2mat(id_dx(:,2));
id_asd = id_all(dx==1);

id_cluster = id_asd(cluster_idx);

%% Get participants' scores
[~, ~, ~, ~, ~, ~, id_scores] = CBIG_ASDf_getSubData(sub_info_file, ...
id_scoreGrp, behav_score_file, idx_score);
scores = cell2mat(id_scores(:,2:end));

%% Get regressors
regressors = CBIG_ASDf_genRegressors(sub_info_file, id_scoreGrp);

%% y
y = ismember(id_scoreGrp, id_cluster); % 1: belong to current cluster; 0: not belong to current cluster

%% Run logistic regression
X = [scores regressors];
n = ones(size(y, 1), 1);

% normalize X
mean_X = mean(X, 1);
std_X = std(X, 0, 1);
X_z = bsxfun(@minus, X, mean_X);
X_z = bsxfun(@rdivide, X_z, std_X);

model = fitglm(X_z, [y n], 'Distribution', 'binomial', 'link', 'logit');
b = table2array(model.Coefficients(:,1)); %coefficients

fprintf('done\n');

disp('Exact LR test');
% Unrestricted model
b_u = b;
yFit = glmval(b_u, X_z, 'logit', 'size', n);
ll_u = sum(log(binopdf(y, n, yFit./n)));

% Overall
% Restricted model
X_r = X_z(:,size(scores,2)+1:end);
%X_r = [];
%dof = size(scores,2);
dof = abs(size(X_z,2) - size(X_r,2)); % degrees of freedom
disp('Overall');
p_lr = CBIG_lr_test(X_r, y, n, ll_u, dof);
fprintf('p = %e\n', p_lr);

%% Remove paths
rmpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','functions'));
rmpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','characteristics'));
rmpath(fullfile(ZHANG_DIR,'step3_analyses_internalUse','MCI2ADProgression'));
