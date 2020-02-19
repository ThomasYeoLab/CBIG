function p_behavior =  CBIG_ASDf_clusterBehavAssoc_wrapper(outputDir)
% p_behavior =  CBIG_ASDf_clusterBehavAssoc_wrapper(outputDir)
% 
% Wrapper function to perform 4 sets of behavioral analyses (only those significant in 
% factor behavioral analyses) on k-means clusters for comparison to behavioral associations in factors. 
%
% Input:
%     - outputDir:
%           Absolute path to the output directory
%
% Output:
%     - p_behavior:
%           All p-values obtained from CCA behavioral association analyses. 
%           These p-values can be used for FDR multiple comparisons correction later.
%
% Example:
%       p_behavior =  CBIG_ASDf_clusterBehavAssoc_wrapper('~/output/kmeans')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
addpath(fullfile(CODE_DIR,'step3_analyses','behavioralAssociation'));
addpath(fullfile(CBIG_CODE_DIR,'external_packages','matlab','non_default_packages','palm','palm-alpha109'));

%% Some constants
kmeans_results = fullfile(UNIT_TEST_DIR,'results','kmeans','kmeans_clusters.mat');
load(kmeans_results); % idx: cluster indices
cluster2factorIdx = [2, 1, 3]; % re-order to match kmeans clusters with factors

input_dir = fullfile(UNIT_TEST_DIR,'data');

p_behavior = [];
%% Prepare inputs
behav_score_file = fullfile(input_dir,'behavior_scores_654.csv');
sub_info_file = fullfile(input_dir,'subInfo_654.csv');

T = readtable(behav_score_file);
score_names = T.Properties.VariableNames;
score_names = CBIG_ASDf_renameBehavScores(score_names);

% Get all subjects' diagnoses
[~, id_dx, ~, ~, ~, ~, ~, ~] = CBIG_ASDf_getSubData(sub_info_file);
id_all = id_dx(:,1);
dx_info = cell2mat(id_dx(:,2));

regressors = CBIG_ASDf_genRegressors(sub_info_file);

%% Restricted/repetitive behaviors
curr_cluster = cluster2factorIdx(1);
fprintf('-----RRB vs cluster %d-----\n', 1);
idx_score = [23:29]; % Indices of behavioral scores are hard-coded
id_scoreGrp = CBIG_ASDf_CCA_genInputs(idx_score,behav_score_file,id_all,dx_info,regressors);

cluster_idx = idx == curr_cluster;
[model, p] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, sub_info_file, behav_score_file, idx_score);

file_name = fullfile(outputDir,'RRB_cluster1');
coef = model.Coefficients(2:length(idx_score)+1,1);
p_scores = model.Coefficients(2:length(idx_score)+1,4);
Rsquared = model.Rsquared;

save(file_name, 'model', 'p', 'coef', 'p_scores', 'Rsquared', 'id_scoreGrp');

%pScores2table(p_scores, score_names(idx_score), coef, file_name);
CBIG_ASDf_CCA_plotBar(score_names(idx_score), table2array(coef), [file_name '_barplot']);

p_behavior = [p_behavior; p];

%% Social responsiveness
curr_cluster = cluster2factorIdx(1);
fprintf('-----Social responsiveness vs cluster %d-----\n', 1);
idx_score = [19:22]; % Indices of behavioral scores are hard-coded
id_scoreGrp = CBIG_ASDf_CCA_genInputs(idx_score,behav_score_file,id_all,dx_info,regressors);

cluster_idx = idx == curr_cluster;
[model, p] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, sub_info_file, behav_score_file, idx_score);

file_name = fullfile(outputDir,'SRS_cluster1');
coef = model.Coefficients(2:length(idx_score)+1,1);
p_scores = model.Coefficients(2:length(idx_score)+1,4);
Rsquared = model.Rsquared;

save(file_name, 'model', 'p', 'coef', 'p_scores', 'Rsquared', 'id_scoreGrp');

%pScores2table(p_scores, score_names(idx_score), coef, file_name);
CBIG_ASDf_CCA_plotBar(score_names(idx_score), table2array(coef), [file_name '_barplot']);

p_behavior = [p_behavior; p];

%%  Comorbid psychopathology (CBCL)
% Cluster 1
curr_cluster = cluster2factorIdx(1);
fprintf('-----Comorbid psychopathology vs cluster %d-----\n', 1);
idx_score = [57:62];
id_scoreGrp = CBIG_ASDf_CCA_genInputs(idx_score,behav_score_file,id_all,dx_info,regressors);

cluster_idx = idx == curr_cluster;
[model, p] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, sub_info_file, behav_score_file, idx_score);

file_name = fullfile(outputDir, 'CBCL_cluster1');
coef = model.Coefficients(2:length(idx_score)+1,1);
p_scores = model.Coefficients(2:length(idx_score)+1,4);
Rsquared = model.Rsquared;

save(file_name, 'model', 'p', 'coef', 'p_scores', 'Rsquared', 'id_scoreGrp');

%pScores2table(p_scores, score_names(idx_score), coef, file_name);
CBIG_ASDf_CCA_plotBar(score_names(idx_score), table2array(coef), [file_name '_barplot']);

p_behavior = [p_behavior; p];

% Cluster 2
curr_cluster = cluster2factorIdx(2);
fprintf('-----Comorbid psychopathology vs cluster %d-----\n', 2);
idx_score = [57:62];
id_scoreGrp = CBIG_ASDf_CCA_genInputs(idx_score,behav_score_file,id_all,dx_info,regressors);

cluster_idx = idx == curr_cluster;
[model, p] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, sub_info_file, behav_score_file, idx_score);

file_name = fullfile(outputDir, 'CBCL_cluster2');
coef = model.Coefficients(2:length(idx_score)+1,1);
p_scores = model.Coefficients(2:length(idx_score)+1,4);
Rsquared = model.Rsquared;

save(file_name, 'model', 'p', 'coef', 'p_scores', 'Rsquared', 'id_scoreGrp');

%pScores2table(p_scores, score_names(idx_score), coef, file_name);
CBIG_ASDf_CCA_plotBar(score_names(idx_score), table2array(coef), [file_name '_barplot']);

p_behavior = [p_behavior; p];

%% Executive function (BRIEF)
curr_cluster = cluster2factorIdx(2);
fprintf('-----Executive function vs cluster %d-----\n', 2);
idx_score = [31:33 35:39];
id_scoreGrp = CBIG_ASDf_CCA_genInputs(idx_score,behav_score_file,id_all,dx_info,regressors);

cluster_idx = idx == curr_cluster;
[model, p] = CBIG_ASDf_logReg_clusterBehavAssoc(cluster_idx, id_scoreGrp, sub_info_file, behav_score_file, idx_score);

file_name = fullfile(outputDir, 'BRIEF_cluster2');
coef = model.Coefficients(2:length(idx_score)+1,1);
p_scores = model.Coefficients(2:length(idx_score)+1,4);
Rsquared = model.Rsquared;

save(file_name, 'model', 'p', 'coef', 'p_scores', 'Rsquared', 'id_scoreGrp');

%pScores2table(p_scores, score_names(idx_score), coef, file_name);
CBIG_ASDf_CCA_plotBar(score_names(idx_score), table2array(coef), [file_name '_barplot']);

p_behavior = [p_behavior; p];

%% Remove path
close all;
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
rmpath(fullfile(CODE_DIR,'step3_analyses','behavioralAssociation'));
rmpath(fullfile(CBIG_CODE_DIR,'external_packages','matlab','non_default_packages','palm','palm-alpha109'));

end
