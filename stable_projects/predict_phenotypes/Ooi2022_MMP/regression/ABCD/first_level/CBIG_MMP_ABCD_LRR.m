function CBIG_MMP_ABCD_LRR(sites, innerFolds, feature_path, featurebase, outdir, ...
    subtxt, subcsv, predvar, covtxt, ymat, covmat, y_idx)

%function CBIG_MMP_ABCD_LRR(sites, innerFolds, feature_path, featurebase, outdir, ...
%   subtxt, subcsv, predvar, covtxt, ymat, covmat, y_idx)
%
% This function prepares the input parameters for the linear ridge
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
%     The full path of output directory where `LRR_(featurebase)` will be created.
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
%   - y_idx
%     The index of the y variable to run the regression for.
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
outstem = convertStringsToChars(strcat("LRR_", featurebase));
feature_name = featurebase;
param.outstem = outstem;
% subject details
sub_txt = subtxt;
sub_csv = {subcsv};
pred_var_txt = predvar;
cov_txt = covtxt;
y_mat = ymat;
cov_mat = covmat;
% regression settings
param.gpso_dir = fullfile(getenv('CBIG_CODE_DIR'), 'external_packages/matlab/non_default_packages/Gaussian_Process');
param.eval = 15;
param.tree = 3;
param.metric = 'predictive_COD';

% modify domain and outstem based on feature
if contains(feature_name,"features_tbss") || contains(feature_name,"features_schaefer") || ...
        strcmp(feature_name,"features_rs") || strcmp(feature_name,"features_mid") || ...
        strcmp(feature_name,"features_sst") || strcmp(feature_name,"features_nback")
    domain = [0.001 0.2; -3 3];
    param.outdir = fullfile(outdir, outstem, strcat('results_var_',num2str(y_idx),...
        '_fs_',num2str(domain(1)),'-',num2str(domain(3)),...
        '_lambda_',num2str(domain(2)),'-',num2str(domain(4))));
else
    domain = [-3 3];
    param.outdir = fullfile(outdir, outstem, strcat('results_var_',num2str(y_idx),...
        '_fs_1',...
        '_lambda_',num2str(domain(1)),'-',num2str(domain(2))));
end
param.domain = domain

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
    cov_types = {'categorical', 'continuous'}; % define covariate types
    cov = CBIG_generate_covariates_from_csv(sub_csv, 'subjectkey', cov_names, cov_types, ...
        sub_txt, 'none', 'none', fullfile(outdir,cov_mat), ',');
else
    fprintf('Using existing covariate file \n')
    cov_temp = load(fullfile(outdir,cov_mat));
    cov = cov_temp.covariates;
end
param.covariates = cov;

%% set other params for setup
fprintf('[4] Loading features... \n')
% load features
features = load(strcat(feature_path ,'.mat'));
features.(feature_name)(features.(feature_name) == 0) = 1e-20;
param.feature_mat = features.(feature_name);

%% LRR workflow
fprintf('[5] Run LRR workflow... \n')
CBIG_LRR_fitrlinear_workflow_1measure(param);

% remove utility path
rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'));
