function p_behavior = CBIG_ASDf_CCA_wrapper(outputDir)
% p_behavior = CBIG_ASDf_CCA_wrapper(outputDir)
% 
% Wrapper function for behavioral analyses in ABIDE-II+GENDAAR datasets using CCA
%
% Input:
%     - outputDir:
%           Absolute path to the output directory
% Output:
%     - p_behavior:
%           All p-values obtained from CCA behavioral association analyses. 
%           These p-values can be used for FDR multiple comparisons correction later.
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

p_behavior = [];

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
addpath(fullfile(CODE_DIR,'step3_analyses','behavioralAssociation'));
addpath(fullfile(CBIG_CODE_DIR,'external_packages','matlab','non_default_packages','palm','palm-alpha109'));

CBIG_REPDATA_DIR = genev('CBIG_REPDATA_DIR');
unit_test_dir=fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
input_dir = fullfile(unit_test_dir,'data');
ref_dir = fullfile(unit_test_dir,'results');

%% Prepare CCA inputs
path_behav_scores = fullfile(input_dir,'behavior_scores_654.csv');
path_subInfo = fullfile(input_dir,'subInfo_654.csv');

T = readtable(path_behav_scores);
score_names = T.Properties.VariableNames;
score_names = CBIG_ASDf_renameBehavScores(score_names);

% Get all subjects' diagnoses
[~, id_dx, ~, ~, ~, ~, ~, ~] = CBIG_ASDf_getSubData(path_subInfo);
id_all = id_dx(:,1);
dx_info = cell2mat(id_dx(:,2));

% Get regressors (age, sex, head motion and sites) of all subjects
[~, regressors_CCA] = CBIG_ASDf_genRegressors(path_subInfo);

Nperm = 10000;

% Set random number seed
rng('default');

%% Restricted/repetitive behaviors
disp('----------Restricted/repetitive behaviors:');
ind_scores = [23:29]; % Indices of behavioral scores are hard-coded
[id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores,path_behav_scores,id_all,dx_info,regressors_CCA);
PAPset=palm_quickperms([],EB,Nperm);

num_factors = 3;
run_num = 94; % Final estimate
factor_order = [3 2 1]; % Reorder the factor
path_factorLoading = fullfile(ref_dir,'visualizeFactors',['k' num2str(num_factors)], ...
    ['r' num2str(run_num)],'factorComp.txt');
%%% CCA for each factor
for factor_idx = 1:num_factors
    num_factors
    factor_idx
    [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior(ind_scores, id_keep, reg, ...
        PAPset, path_subInfo, path_behav_scores, path_factorLoading, factor_order, factor_idx);
    % Plot the results
    if Ncca > 0
        file_name = fullfile(outputDir, ['RRB_k=3_F' num2str(factor_idx)]);
        CBIG_ASDf_CCA_plotBar(score_names(ind_scores), strucCorr_score, [file_name '_barPlot']);
        CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx, [file_name '_scatterPlot']);
    end
    p_behavior = [p_behavior; pVal];
end

%% Social responsiveness (SRS)
disp('----------Social responsiveness:');
ind_scores = [19:22]; %hard-coded
[id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores,path_behav_scores,id_all,dx_info,regressors_CCA);
PAPset=palm_quickperms([],EB,Nperm);

num_factors = 3;
run_num = 94; % Final estimate
factor_order = [3 2 1]; % Reorder the factor
path_factorLoading = fullfile(ref_dir,'visualizeFactors',['k' num2str(num_factors)], ...
    ['r' num2str(run_num)],'factorComp.txt');
%%% CCA for each factor
for factor_idx = 1:num_factors
    num_factors
    factor_idx
    [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior(ind_scores, id_keep, reg, ...
        PAPset, path_subInfo, path_behav_scores, path_factorLoading, factor_order, factor_idx);
    % Plot the results
    if Ncca > 0
        file_name = fullfile(outputDir, ['SRS_k=3_F' num2str(factor_idx)]);
        CBIG_ASDf_CCA_plotBar(score_names(ind_scores), strucCorr_score, [file_name '_barPlot']);
        CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx, [file_name '_scatterPlot']);
    end
    p_behavior = [p_behavior; pVal];
end


%% Comorbid psychopathology (CBCL)
disp('----------Comorbid psychopathology:');
ind_scores = [57:62];
[id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores,path_behav_scores,id_all,dx_info,regressors_CCA);
PAPset=palm_quickperms([],EB,Nperm);

num_factors = 3;
run_num = 94; % Final estimate
factor_order = [3 2 1]; % Reorder the factor
path_factorLoading = fullfile(ref_dir,'visualizeFactors',['k' num2str(num_factors)], ...
    ['r' num2str(run_num)],'factorComp.txt');
%%% CCA for each factor
for factor_idx = 1:num_factors
    num_factors
    factor_idx
    [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior(ind_scores, id_keep, reg, ...
        PAPset, path_subInfo, path_behav_scores, path_factorLoading, factor_order, factor_idx);
    % Plot the results
    if Ncca > 0
        file_name = fullfile(outputDir, ['CBCL_k=3_F' num2str(factor_idx)]);
        CBIG_ASDf_CCA_plotBar(score_names(ind_scores), strucCorr_score, [file_name '_barPlot']);
        CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx, [file_name '_scatterPlot']);
    end
    p_behavior = [p_behavior; pVal];
end

%% Executive function (BRIEF)
disp('----------Executive function:');
ind_scores = [31:33 35:39];
[id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores,path_behav_scores,id_all,dx_info,regressors_CCA);
PAPset=palm_quickperms([],EB,Nperm);

%%% K=3
num_factors = 3;
run_num = 94; % Final estimate
factor_order = [3 2 1]; % Reorder the factor
path_factorLoading = fullfile(ref_dir,'visualizeFactors',['k' num2str(num_factors)], ...
    ['r' num2str(run_num)],'factorComp.txt');
%%% CCA for each factor
for factor_idx = 1:num_factors
    num_factors
    factor_idx
    [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior(ind_scores, id_keep, reg, ...
        PAPset, path_subInfo, path_behav_scores, path_factorLoading, factor_order, factor_idx);
    % Plot the results
    if Ncca > 0
        file_name = fullfile(outputDir, ['BRIEF_k=3_F' num2str(factor_idx)]);
        CBIG_ASDf_CCA_plotBar(score_names(ind_scores), strucCorr_score, [file_name '_barPlot']);
        CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx, [file_name '_scatterPlot']);
    end
    p_behavior = [p_behavior; pVal];
end

%% Autistic symptoms (ADOS)
disp('-----------Autistic symptoms:');
ind_scores = [11:13];
[id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores,path_behav_scores,id_all,dx_info,regressors_CCA);
PAPset=palm_quickperms([],EB,Nperm);

%%% K=3
num_factors = 3;
run_num = 94; % Final estimate
factor_order = [3 2 1]; % Reorder the factor
path_factorLoading = fullfile(ref_dir,'visualizeFactors',['k' num2str(num_factors)], ...
    ['r' num2str(run_num)],'factorComp.txt');
%%% CCA for each factor
for factor_idx = 1:num_factors
    num_factors
    factor_idx
    [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior(ind_scores, id_keep, reg, ...
        PAPset, path_subInfo, path_behav_scores, path_factorLoading, factor_order, factor_idx);
    % Plot the results
    if Ncca > 0
        file_name = fullfile(outputDir, ['ADOS_k=3_F' num2str(factor_idx)]);
        CBIG_ASDf_CCA_plotBar(score_names(ind_scores), strucCorr_score, [file_name '_barPlot']);
        CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx, [file_name '_scatterPlot']);
    end
    p_behavior = [p_behavior; pVal];
end

%% Remove path
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
rmpath(fullfile(CODE_DIR,'step3_analyses','behavioralAssociation'));
rmpath(fullfile(CBIG_CODE_DIR,'external_packages','matlab','non_default_packages','palm','palm-alpha109'));

end
