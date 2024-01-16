function CBIG_MMP_HCP_Elasticnet(outerFolds, innerFolds, feature_path, featurebase, split, outdir, ...
    subtxt, scorecsv, restrictedcsv, predvar, covtxt, ymat, covmat, y_idx)

% function CBIG_MMP_HCP_Elasticnet(outerFolds, innerFolds, feature_path, featurebase, split, outdir, ...
%    subtxt, scorecsv, restrictedcsv, predvar, covtxt, ymat, covmat, y_idx)
%
% This function prepares the input parameters for a specified cross-validated split
% of elasticnet regression.
%
% Inputs:
%   - outerFolds
%     Number of outer-loop cross validation folds.
%
%   - innerFolds
%     Number of inner-loop cross validation folds.
%
%   - feature_path
%     Full path to the feature file used for the regression. This function assumes that
%     feature matrices are saved in a mat file with a struct containing a
%     #features x # subjects 2-D matrix. The matrix should be named as `featurebase` (see below).
%
%   - featurebase
%     Base file name for the input features. Output folders generated in a folder `KRR_(featurebase)`.
%
%   - split
%     Random seed for split to be generated.
%
%   - outdir
%     The full path of output directory where `KRR_(featurebase)` will be created.
%
%   - subtxt
%     Full path of the subject ID list. Each line in this list corresponds
%     to one subject ID. This list should be in the same order as the feature
%     matrix.
%
%   - scorecsv
%     Full path of csv file containing behavioral data and gender for all subjects.
%
%   - restrictedcsv
%     Full path of csv file containing restricted data (Family ID and age) for all subjects.
%
%   - predvar
%     Full path to a text file with all behavioral (or demographic)
%     measures (measures to be predicted using kernel ridge regression).
%     Each line corresponds to one behavioral name. The behavioral names
%     should exist as a header in "scorecsv".
%
%   - covtxt
%     Full path to a text file stating all covariate names. Each line
%     corresponds to one covariate name. The covariate names should exist
%     as header in the "restrictedcsv".
%
%   - ymat
%     Name of the output mat file containing behavioural variables for all
%     subjects in the regression analysis.
%
%   - covmat
%     Name of the output mat file containing covariates for all subjects in the
%     regression analysis.
%
%   - y_idx
%     Index of y variable.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% format params that are read in
% folds
seed = split;
seed_name = strcat('seed_', num2str(seed));
num_outer_folds = outerFolds;
param.num_inner_folds = innerFolds;
% feature details
outstem = convertStringsToChars(strcat("Elasticnet_", featurebase));
feature_name = featurebase;
param.outstem = outstem;
param.split_name = seed_name;
param.outdir = fullfile(outdir, outstem, strcat('results_var_',num2str(y_idx)));
% subject details
sub_txt = subtxt;
score_csv = {scorecsv}; % eg. diff_ct_intersection.csv
restricted_csv = {restrictedcsv};
cov_csv = {scorecsv restrictedcsv}; % eg. diff_ct_intersection.csv
pred_var_txt = predvar;
cov_txt = covtxt;
y_mat = ymat;
cov_mat = covmat;
% regression settings
lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
alpha_set = [0.01 0.1 0.4 0.7 0.9 0.99];
param.lambda = lambda_set;
param.alpha = alpha_set;
param.metric = 'predictive_COD';

%% get subfold
fprintf('[1] Generate subfold... \n')

% generate folds
fold_mat = 'no_relative_10_fold_sub_list.mat';
if ~exist(fullfile(param.outdir, fold_mat))
    sub_fold = CBIG_cross_validation_data_split(sub_txt, restricted_csv{1}, ...
        'Subject', 'Family_ID', num_outer_folds, seed, param.outdir, ',');
else
    fprintf('Using existing sub_fold file \n')
    fold_temp = load(fullfile(param.outdir, fold_mat));
    sub_fold = fold_temp.sub_fold;
end
param.sub_fold  = sub_fold;

%% generate y matrix
fprintf('[2] Generate y matrix... \n')

if ~exist(fullfile(outdir, y_mat))
    % get names of tasks to predict
    fid = fopen(pred_var_txt,'r'); % variable names text file
    score_list = textscan(fid,'%s');
    score_names = score_list{1};
    fclose(fid);
    num_scores = size(score_names,1);
    % generate y
    score_types = cell(1,num_scores); % define score types
    score_types(:) = {'continuous'};
    y = CBIG_read_y_from_csv(score_csv, 'Subject', score_names, score_types,...
        sub_txt, fullfile(outdir, y_mat), ',');
else
    fprintf('Using existing y file \n')
    y_temp = load(fullfile(outdir,y_mat));
    y = y_temp.y;
end
param.y = y;

% choose y to run
param.y = y(:,y_idx);

%% generate covariate matrix
fprintf('[3] Generate covariate matrix... \n')

if ~exist(fullfile(outdir,cov_mat))
    % generate covariates
    fid = fopen(cov_txt,'r'); % covariate names text file
    cov_list = textscan(fid,'%s');
    cov_names = cov_list{1};
    fclose(fid);
    num_cov = size(cov_names,1);
    cov_types = {'continuous', 'categorical'}; % define covariate types
    cov = CBIG_generate_covariates_from_csv(cov_csv, 'Subject', cov_names, cov_types, ...
        sub_txt, 'none', 'none', fullfile(outdir,cov_mat), ',');
else
    fprintf('Using existing covariate file \n')
    cov_temp = load(fullfile(outdir, cov_mat));
    cov = cov_temp.covariates;
end
param.covariates = cov;

%% set other params for setup
fprintf('[4] Loading features... \n')
% load features
features = load(strcat(feature_path ,'.mat'));
param.feature_mat = features.(feature_name);

%% Elasticnet workflow
fprintf('[5] Run Elasticnet workflow... \n')
CBIG_run_Elasticnet_workflow(param);
