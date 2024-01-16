function CBIG_MMP_ABCD_KRR(sites, innerFolds, feature_path, featurebase, outdir, ...
    subtxt, subcsv, predvar, covtxt, ymat, covmat)

%function CBIG_MMP_ABCD_KRR(sites, innerFolds, feature_path, featurebase, outdir, ...
%   subtxt, subcsv, predvar, covtxt, ymat, covmat)
%
% This function prepares the input parameters for the single-kernel
% regression leave-p-out workflow.
%
% Inputs:
%   - sites
%     Number of sites used in the test fold for the leave-p-out cross validation.
%
%   - innerFolds
%     Number of inner-loop cross validation folds
%
%   - feature_path
%     Full path to the feature file used for the regression. This function assumes that
%     feature matrices are saved in a mat file with a struct containing a
%     #features x # subjects 2-D matrix. The matrix should be named as `featurebase` (see below).
%
%   - featurebase
%     Base file name for the input features. Output folders generated in a folder `KRR_(featurebase)`.
%
%   - outdir
%     The full path of output directory where `KRR_(featurebase)` will be created.
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
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add utility path
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'));

%% format params that are read in
% folds
num_sites = sites;
param.num_inner_folds = innerFolds;
% feature details
outstem = convertStringsToChars(strcat("KRR_", featurebase));
feature_name = featurebase;
param.outstem = outstem;
param.outdir = fullfile(outdir, outstem, 'results');
% subject details
sub_txt = subtxt;
sub_csv = {subcsv};
pred_var_txt = predvar;
cov_txt = covtxt;
y_mat = ymat;
cov_mat = covmat;
% regression settings
lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
param.with_bias = 1;
param.ker_param.type = 'corr';
param.ker_param.scale = NaN;
param.lambda_set = lambda_set;
param.threshold_set = [];
param.cov_X = [];
param.metric = 'predictive_COD';

%% get subfold
fprintf('[1] Generate subfold... \n')

% generate folds
fold_mat = 'no_relative_3_fold_sub_list.mat';
if ~exist(fullfile(outdir, fold_mat))
    sub_fold = CBIG_MMP_ABCD_LpOCV_split( sub_txt, sub_csv{1}, ...
        'subjectkey', 'site_group', num_sites, outdir, ',' );
else
    fprintf('Using existing sub_fold file \n')
    fold_temp = load(fullfile(outdir,fold_mat));
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
    y = CBIG_read_y_from_csv(sub_csv, 'subjectkey', score_names, score_types,...
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
    cov_types = {'categorical', 'continuous'}; % define covariate types
    cov = CBIG_generate_covariates_from_csv(sub_csv, 'subjectkey', cov_names, cov_types, ...
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

%% KRR workflow
fprintf('[5] Run KRR workflow...')
CBIG_KRR_workflow(param);

% remove utility path
rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'));