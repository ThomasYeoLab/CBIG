function CBIG_TRBPC_LRR_LpOCV_prepare_parameters( csv_file, subject_list, feature_file, ...
    y_list, covariate_list, FD_file, DVARS_file, outdir, outstem, num_leave_out, num_inner_folds, ...
    domain, eval, tree, metric )

% CBIG_TRBPC_LRR_LpOCV_prepare_parameters( csv_file, subject_list, feature_file, ...
%    y_list, covariate_list, FD_file, DVARS_file, outdir, outstem, num_leave_out, num_inner_folds, ...
%    domain, eval, tree, metric )
%
% This function prepares the input parameters for the linear ridge regression leave-p-out workflow
%
% Inputs:
%   - csv_file
%     Full path of csv file contain behavioral data and covariates for all subejcts.
%
%   - subject_list
%     Full path of the subject ID list. Each line in this list corresponds
%     to one subject ID.
%
%   - feature_files
%     Full path of the file storing a a cell array. The cell array contains
%     the full path of all the feature files that are required to
%     calculate the kernels. For each feature file, a matrix "feature_mat" 
%     is assumed to be saved in this file. "feature_mat" can be a 2-D matrix
%      with dimension of #features x #subjects or a 3-D
%     matrix with dimension of #ROIs1 x #ROIs2 x #subjects. If
%     "feature_mat" is 3-D, it is a connectivity matrix between two sets of
%     ROIs, and only the lower-triangular off-diagonal entries will be
%     considered as features because the connectivity matrix is symmetric.
%
%   - y_list
%     Full path to a text file with all behavioral (or demographic)
%     measures (measures to be predicted using kernel ridge regression).
%     Each line corresponds to one behavioral name. The behavioral names
%     should exist as a header in csv_file
%
%   - covariate_list
%     Full path to a text file stating all covariate names. Each line
%     corresponds to one covariate name. The covariate names should exist
%     as header in cvs_files.
%
%   - FD_file (optional)
%     If there is a need to regress 'FD' from the behavioral (or demographic)
%     measures, y, the user should include 'FD' in the "covariate_list". In
%     this case, "FD_file" is the full path of the mean framewise
%     displacement (FD) of all subjects. The number of lines in "FD_file"
%     should be the same as the number of lines in "subject_list".
%     If the user does not need to regress 'FD' from y, then the input
%     variable 'FD_file' is not required and the user can pass in 'NONE' to
%     the function.
%     If "covariate_list" does not contain FD, this argument will be
%     ignored.
%
%   - DVARS_file (optional)
%     If there is a need to regress 'DVARS' from the behavioral
%     (or demographic) measures, y, the user must include the covariate
%     'DVARS' (or 'DV') in the 'covariate_list'. In this case, "DVARS_file"
%     is the full path of the mean DVARS of all subjects. The number of
%     lines in "DVARS_file" should be the same as the number of lines in
%     "subject_list".
%     If the user does not need to regress 'DVARS' from y, then the input
%     variable 'DVARS_file' is not required and the user can pass in 'NONE'
%     to the function.
%     If "covariate_list" does not contain DV (or DVARS), this argument
%     will be ignored.
%
%   - outdir
%     The full path of output directory.
%
%   - outstem
%     A string appended to the file names to specify the output files
%     (y after regression, accuracy files, ...). For example, if outstem =
%     'allscore', then the accuracy files will be names as
%     <path_to_file>/acc_allscore.mat, and the final output filename will be
%     [outdir '/final_result_allscore.mat'].
%     If no outstem is required, the user can just pass in an empty string
%     ('').
%
%   - num_leave_out
%     A string or scalar, the value of p for the leave-p-out cross-validation.
%
%   - num_inner_folds
%     A string or scalar.
%     To select optimal hyperparameters, each training fold will be split
%     randomly into "num_inner_folds" inner-loop cross-validation folds.
%
%   - domain (optional)
%     Full path of the file (.mat) storing the domain used to 
%     specify the search boundary of the hyperparameter. A 1 x 2 vector
%     "domain" is assumed to be saved in this file.If this file is not
%     passed in and also does NOT exist in
%     "setup_file", or "domain" is 'NONE', it will be set as default: 
%     [0 20]. The first element of the domain stores the lower bound and
%     the second element stores the upper bound
%
%   - metric (optional)
%     A string stating which accuracy metric to be used for optimising
%     hyperparameters. Currently, the following are accepted:
%     Choose from:
%       'corr'              - Pearson's correlation;
%       'COD'               - Coefficient of determination. Defined as
%                             1-||y_pred-y_test||^2/||mean(y_test)-y_test||^2,
%                             where y_pred is the prediction of the test data, 
%                             y_test is the groud truth of the test data, 
%                             and mean(y_test) is the mean of test data
%       'predictive_COD'    - Predictive coefficient of determination. Defined as
%                             1-||y_pred-y_test||^2/||mean(y_train)-y_test||^2,
%                             where y_pred is the prediction of the test data, 
%                             y_test is the groud truth of the test data, 
%                             and mean(y_train) is the mean of training data
%       'MAE'               - mean absolute error
%       'MAE_norm'          - mean absolute error divided by the standard
%                             derivation of the target variable of the training set
%       'MSE'               - mean squared error
%       'MSE_norm'          - mean squared error divided by the variance
%                             of the target variable of the training set
%     If not provided, 'predictive_COD' is default.
%
% Outputs:
%   This function with save out necessary files in the outdir for the multiKRR LpOCV
%   workflow. The files includes:
%   - regression model hyperparameters
%   - prediction target variable y and regressed y
%   - the covariates
%   - cross-validation split indices
%   - kernels for each input feature set
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% setting up
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','predict_phenotypes', 'ChenTam2022_TRBPC');
addpath(genpath(project_code_dir));

if(ischar(num_leave_out))
    num_leave_out = str2double(num_leave_out);
end

if(ischar(num_inner_folds))
    num_inner_folds = str2double(num_inner_folds);
end

%% Data split
fprintf('[LRR LpOCV workflow]: leave-%d-sites-out cross-validation splits.\n',num_leave_out);
CBIG_TRBPC_LpOCV_split( subject_list, csv_file, 'subjectkey', ...
    'site_group_150', num_leave_out, outdir, ',' );

%% Read y
% y types
fprintf('[LRR LpOCV workflow]: read the measures to be predicted.\n')
[y_names, num_y] = CBIG_text2cell(y_list);
for i = 1:num_y
    if(strcmp(y_names{i}, 'gender'))
        y_types{i} = 'categorical';
    else
        y_types{i} = 'continuous';
    end
end
ystem = outstem;
if(~isempty(outstem))
    ystem = ['_' outstem];
end
if(~exist([outdir '/y' ystem '.mat'], 'file'))
    CBIG_read_y_from_csv( {csv_file}, 'subjectkey', y_names, y_types, ...
        subject_list, fullfile(outdir, ['y' ystem '.mat']), ',' );
end

%% Read covariates
% covariate types
fprintf('[LRR LpOCV workflow]: read covariates to be regressed from the measures.\n')
[cov_names, num_cov] = CBIG_text2cell(covariate_list);
for i = 1:num_cov
    if(strcmp(cov_names{i}, 'gender') || strcmp(cov_names{i}, 'race_ethnicity'))
        cov_types{i} = 'categorical';
    else
        cov_types{i} = 'continuous';
    end
end
cov_stem = outstem;
if(~isempty(outstem))
    cov_stem = ['_' outstem];
end
if(~exist(fullfile(outdir, ['covariates' cov_stem '.mat']), 'file'))
    CBIG_generate_covariates_from_csv( {csv_file}, ...
        'subjectkey', cov_names, cov_types, subject_list, FD_file, DVARS_file, ...
        fullfile(outdir, ['covariates' cov_stem '.mat']), ',' );
end

%% generate parameters
fprintf('[LRR LpOCV workflow]: generate input parameters ...\n')
sub_fold_file = fullfile(outdir, ['no_relative_' num2str(num_leave_out) '_fold_sub_list.mat']);
param = CBIG_LRRworkflow_parseInput( sub_fold_file, fullfile(outdir, ['y' ystem '.mat']), ...
    fullfile(outdir, ['covariates' cov_stem '.mat']), num_inner_folds, ...
    outdir, outstem, domain, eval, tree, metric);
% load feature matrix
feature_struct = load(feature_file);
feature_name = fieldnames(feature_struct);
param.FC_mat = feature_struct.(feature_name{1});

save(fullfile(outdir, 'param.mat'),'param','-v7.3');

rmpath(genpath(project_code_dir));
end

function param = CBIG_LRRworkflow_parseInput(varargin)

varnames = {'sub_fold', 'y', 'covariates', 'num_inner_folds', 'outdir', ...
    'outstem', 'domain', 'eval', 'tree', 'metric'};

for i = 1:(numel(varnames)-4)
    if(length(varargin)<i || isempty(varargin{i}))
        error('Variable ''%s'' is needed.', varnames{i})
    else
        if(i~=4 && i~=5 && i~=6)
            % sub_fold, y, covariates, feature_mat are loaded from .mat files.
            currvar = load(varargin{i});
            names = fieldnames(currvar);
            if(any(strcmp(names, varnames{i})))
                currname = varnames{i};
            else
                currname = names{1};
            end
            param.(varnames{i}) = currvar.(currname);
        else
            % num_inner_folds or outdir are directly set as fields of param
            param.(varnames{i}) = varargin{i};
        end
    end
end

% default of domain
if(length(varargin)<7 || isempty(varargin{7}) || strcmpi(varargin{7}, 'none'))
    param.domain = [0.001 0.1; 3 8];
else
    load(varargin{7})
    param.domain = domain;
end

% default of eval
if(length(varargin)<8 || isempty(varargin{8}) || strcmpi(varargin{8}, 'none'))
    param.eval = 15;
else
    param.domain = varargin{8};
end

% default of tree
if(length(varargin)<9 || isempty(varargin{9}) || strcmpi(varargin{9}, 'none'))
    param.tree = 3;
else
    param.tree = varargin{9};
end

% default of metric
if(length(varargin)<10 || isempty(varargin{10}))
    param.metric = 'predictive_COD';
else
    param.metric = varargin{10};
end

end
