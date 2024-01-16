function CBIG_MMP_HCP_multiKRR_example(outerFolds, innerFolds, feature_dir, feature_mats_cell, ...
    kernel_groups, outstem_name, split_idx, outdir, subtxt, scorecsv, restrictedcsv, ...
    predvar, covtxt, ymat, covmat, fold_idx)

% function CBIG_MMP_HCP_multiKRR(outerFolds, innerFolds, feature_dir, feature_mats_cell, ...
%    kernel_groups, outstem_name, split_idx, outdir, subtxt, subcsv, predvar, ...
%    covtxt, ymat, covmat, fold_idx)
%
% This function prepares the input parameters for the multi-kernel
% regression cross-validation workflow.
%
% Inputs:
%   - sites
%     Number of sites used in the test fold for the leave-p-out cross validation.
%
%   - innerFolds
%     Number of inner-loop cross validation folds
%
%   - feature_dir
%     Path to a folder containing all feature files used for the regression. This function assumes
%     that feature matrices are saved in a mat file with a struct containing a
%     #features x # subjects 2-D matrix. The matrix should be named as `featurebase` (see below).
%
%   - feature_mats_cell
%     Cell of base file name for the input features.
%
%   - kernel_groups
%     Kernel groupings for multiKRR workflow. For example if 4 kernels are to be used to train the
%     multi kernel model, but they can be grouped such that 3 of the kernels uses the same hyperparameter
%     lambda for regularization, then I can pass in the 'group_kernel' as {[1,3,4],[2]}. This means that
%     kernels 1, 3 and 4 are regularized using the same lambda while kernel 2 is regularized using a different
%     lambda. Do note that kernel 1 would correspond to the first feature mat in 'feature_file' so on and so forth.
%
%   - outstem_name
%     Name of output folder.
%
%   - split_idx
%     The split which the multiKRR workflow is being run for.
%
%   - outdir
%     The full path of output directory where `outstem_name` will be created.
%
%   - subtxt
%     Full path of the subject ID list. Each line in this list corresponds
%     to one subject ID. This list should be in the same order as the feature
%     matrix.
%
%   - subcsv
%     Full path of csv file containing behavioral data and covariates for all subjects.
%
%   - predvar
%     Full path to a text file with all behavioral (or demographic)
%     measures (measures to be predicted using kernel ridge regression).
%     Each line corresponds to one behavioral name. The behavioral names
%     should exist as a header in "subcsv".
%
%   - covtxt
%     Full path to a text file stating all covariate names. Each line
%     corresponds to one covariate name. The covariate names should exist
%     as header in the "subcsv".
%
%   - ymat
%     Name of the output mat file containing behavioural variables for all
%     subjects in the regression analysis.
%
%   - covmat
%     Name of the output mat file containing covariates for all subjects in the
%     regression analysis.
%
%   - fold_idx
%     The fold which the multiKRR workflow is being run for.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add utility path
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'utilities'));

%% format params that are read in
% folds
seed = split_idx;
seed_name = strcat('seed_', num2str(seed));
num_outer_folds = outerFolds;
param.num_inner_folds = innerFolds;
% feature details
outstem = outstem_name;
param.outstem = outstem;
param.outdir = fullfile(outdir, outstem, seed_name, 'results')
% subject details
sub_txt = subtxt;
score_csv = {scorecsv};
%restricted_csv = {restrictedcsv};
%cov_csv = {scorecsv restrictedcsv};
pred_var_txt = predvar;
cov_txt = covtxt;
y_mat = ymat;
cov_mat = covmat;
% regression settings
lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
param.lambda_set = lambda_set;
param.domain = [0 20];
param.with_bias = 1;
param.ker_param.type = 'corr';
param.ker_param.scale = NaN;
param.threshold = 'None';
param.group_kernel = kernel_groups;
param.acc_metric = 'predictive_COD';

%% get subfold
fprintf('[1] Generate subfold... \n')

% generate folds
fold_mat = 'no_relative_3_fold_sub_list.mat';
if ~exist(fullfile(param.outdir, fold_mat))
    sub_fold = CBIG_cross_validation_data_split(sub_txt, 'NONE', ...
        'NONE', 'NONE', num_outer_folds, seed, param.outdir, ',');
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

%% generate covariate matrix
fprintf('[3] Generate covariate matrix... \n')

if ~exist(fullfile(outdir,cov_mat))
    % generate covariates
    fid = fopen(cov_txt,'r'); % covariate names text file
    cov_list = textscan(fid,'%s');
    cov_names = cov_list{1};
    fclose(fid);
    num_cov = size(cov_names,1);
    cov_types = {'categorical'}; % define covariate types
    cov = CBIG_generate_covariates_from_csv(score_csv, 'Subject', cov_names, cov_types, ...
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
for n = 1:size(feature_mats_cell,2)
    feature_name = feature_mats_cell{n};
    % create cell with full paths
    feature_paths.feature_mat{n} = fullfile(feature_dir, strcat(feature_name ,'.mat'));
end
feat_mat_file = fullfile(outdir, outstem, 'feat_mat_path.mat');
save(feat_mat_file, '-struct', 'feature_paths');

%% Multi-KRR workflow
fprintf('[5] Run multi-KRR workflow...\n')
CBIG_MMP_MultiKRR_workflow_parallel_subfolds('',0, fullfile(param.outdir,fold_mat), ...
    fullfile(outdir,y_mat), fullfile(outdir,cov_mat), feat_mat_file, ...
    param.num_inner_folds, param.outdir, param.outstem, fold_idx);

% remove utility path
rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'utilities'));
