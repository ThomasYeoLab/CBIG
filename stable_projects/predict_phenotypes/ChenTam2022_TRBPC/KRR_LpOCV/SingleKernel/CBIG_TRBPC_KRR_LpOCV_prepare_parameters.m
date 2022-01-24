function CBIG_TRBPC_KRR_LpOCV_prepare_parameters( csv_file, subject_list, feature_file, ...
    y_list, covariate_list, FD_file, DVARS_file, outdir, outstem, num_leave_out, num_inner_folds, ...
    ker_param_file, lambda_set_file, threshold_set_file, metric )

%function CBIG_TRBPC_KRR_LpOCV_prepare_parameters( csv_file, subject_list, feature_file, ...
%    y_list, covariate_list, FD_file, DVARS_file, outdir, outstem, num_leave_out, num_inner_folds, ...
%    ker_param_file, lambda_set_file, threshold_set_file, metric )
% 
% This function prepares the input parameters for the single-kernel
% regression leave-p-out workflow
% 
% Inputs:
%   - csv_file
%     Full path of csv file contain behavioral data and covariates for all subejcts.
%
%   - subject_list
%     Full path of the subject ID list. Each line in this list corresponds
%     to one subject ID.
%  
%   - feature_file
%     Full path of the feature file used to calculate the kernels. A matrix
%     "feature_mat" is assumed to be saved in this file. "feature_mat" can
%     be a 2-D matrix with dimension of #features x #subjects or a 3-D
%     matrix with dimension of #ROIs1 x #ROIs2 x #subjects. If
%     "feature_mat" is 3-D, it is a connectivity matrix between two sets of
%     ROIs, and only the lower-triangular off-diagonal entries will be
%     considered as features because the connectivity matrix is symmetric.
%     If this file is not passed in, then the matrix "feature_mat" is
%     assumed to be stored in "setup_file".
%  
%   - y_list
%     Full path to a text file with all behavioral (or demographic)
%     measures (measures to be predicted using kernel ridge regression). 
%     Each line corresponds to one behavioral name. The behavioral names
%     should exist as a header in the "csv_file"
% 
%   - covariate_list
%     Full path to a text file stating all covariate names. Each line
%     corresponds to one covariate name. The covariate names should exist
%     as header in the "csv_file", except for 'FD' and 'DVARS'.
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
%     The full path of output directory. A subfolder 
%     [outdir '/randseed_' seed] will be created to save all output files 
%     generated using the current random seed.
% 
%   - outstem
%     A string appended to the file names to specify the output files 
%     (y after regression, accuracy files, ...). For example, if outstem =
%     '58behaviors', then the accuracy files will be names as
%     <path_to_file>/acc_58behaviors.mat,
%     and the final output filename will be 
%     [outdir '/randseed_' seed '/final_result_58behaviors.mat'].
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
% 
%   - ker_param_file (optional)
%     Full path of the kernel parameter file (.mat). A structure "ker_param" 
%     is assumed to be saved in this file.
%     "ker_param" is a K x 1 structure with two fields: type and scale. K
%     denotes the number of kernels.
%     ker_param(k).type is a string of the type of k-th kernel. Choose from
%                       'corr'        - Pearson's correlation;
%                       'Gaussian'    - Gaussian kernel;
%                       'Exponential' - exponential kernel.
%     ker_param(k).scale is a scalar specifying the scale of k-th kernel
%     (for Gaussian kernel or exponential kernel). If ker_param(k).type == 'corr',
%     ker_param(k).scale = NaN.
%     If this argument is not passed in (or passed in as 'NONE'), then
%     ker_param will be set as default:
%     ker_param.type = 'corr';
%     ker_param.scale = NaN.
%  
%   - lambda_set_file (optional)
%     Full path of the regularization parameter file (.mat). A vector 
%     "lambda_set" is assumed to be saved in this file.
%     "lambda_set" is a vector of numbers for grid search of lambda (the
%     regularization parameter). If this file is not passed in (or passed
%     in as 'NONE'), it will be set as default:
%     [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
%        5 10 15 20]
% 
%   - threshold_set_file (optional)
%     Full path of the file (.mat) storing the set of threshold used to 
%     binarize the predicted score when the original y is binary. A vector
%     "threshold_set" is assumed to be saved in this file.
%     "threshold_set" is a vector used for grid search of optimal
%     "threshold". If this file is not passed in (or passed in as 'NONE'),
%     or "threshold_set" is 'NONE', it will be set as default:
%     [-1:0.1:1].
%   
%   - metric
%     A string indicating the metric used to define prediction loss. The
%     loss is used to choose hyperparameters.
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

if(~exist('ker_param_file', 'var') || isempty(ker_param_file) || ...
        strcmpi(ker_param_file, 'none'))
    ker_param_file = [];
end

if(~exist('lambda_set_file', 'var') || isempty(lambda_set_file) || ...
        strcmpi(lambda_set_file, 'none'))
    lambda_set_file = [];
end

if(~exist('threshold_set_file', 'var') || isempty(threshold_set_file) || ...
        strcmpi(threshold_set_file, 'none'))
    threshold_set_file = [];
end

%% Data split
fprintf('[KRR LpOCV workflow]: leave-%d-sites-out cross-validation splits.\n',num_leave_out);
CBIG_TRBPC_LpOCV_split( subject_list, csv_file, 'subjectkey', ...
    'site_group_150', num_leave_out, outdir, ',' );

%% Read y
% y types
fprintf('[KRR LpOCV workflow]: read the measures to be predicted.\n')
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
fprintf('[KRR LpOCV workflow]: read covariates to be regressed from the measures.\n')
cov_stem = outstem;
if(~isempty(outstem))
    cov_stem = ['_' outstem];
end

[cov_names, num_cov] = CBIG_text2cell(covariate_list);
if length(cov_names)==1 && strcmpi(cov_names{1},'none')
    covariates = [];
    save(fullfile(outdir, ['covariates' cov_stem '.mat']),'covariates')
else
    for i = 1:num_cov
        if(strcmp(cov_names{i}, 'gender') || strcmp(cov_names{i}, 'race_ethnicity'))
            cov_types{i} = 'categorical';
        else
            cov_types{i} = 'continuous';
        end
    end
    if(~exist(fullfile(outdir, ['covariates' cov_stem '.mat']), 'file'))
        CBIG_generate_covariates_from_csv( {csv_file}, ...
            'subjectkey', cov_names, cov_types, subject_list, FD_file, DVARS_file, ...
            fullfile(outdir, ['covariates' cov_stem '.mat']), ',' );
    end
end

%% generate parameters
sub_fold_file = fullfile(outdir, ['no_relative_' num2str(num_leave_out) '_fold_sub_list.mat']);
param = CBIG_TRBPC_KRR_LpOCV_parseInput( sub_fold_file, fullfile(outdir, ['y' ystem '.mat']), ...
    fullfile(outdir, ['covariates' cov_stem '.mat']), num_inner_folds, ...
    outdir, outstem,1, ker_param_file, lambda_set_file, threshold_set_file,metric);

if(ischar(param.with_bias))
    param.with_bias = str2num(param.with_bias);
end

for col = 1:size(param.y, 2)
    bin_flag(col) = numel(unique(param.y(:,col)))==2;
end
if(any(bin_flag==1) && any(bin_flag==0))
    error('Mixture of binary (e.g. sex) and continuous cases. Please run them separately.')
elseif(any(bin_flag==1))
    param.bin_flag = 1;
else
    param.bin_flag = 0;
    param.threshold_set = NaN;
end

save(fullfile(outdir, 'param.mat'),'param');

%% regress covariates from y
fprintf('[KRR LpOCV workflow]: regress covariates from y for each fold.\n')
CBIG_crossvalid_regress_covariates_from_y( ...
    param.y, param.covariates, param.sub_fold, param.outdir, param.outstem);

%% generate kernels
fprintf('[KRR LpOCV workflow]: generate kernels.\n')
feature_struct = load(feature_file);
feature_name = fieldnames(feature_struct);
feature_mat = feature_struct.(feature_name{1});
CBIG_KRR_generate_kernels_LITE(feature_mat, param.outdir,1, param.ker_param )

rmpath(genpath(project_code_dir));
end

function param = CBIG_TRBPC_KRR_LpOCV_parseInput(varargin)

varnames = {'sub_fold', 'y', 'covariates', 'num_inner_folds', 'outdir', ...
    'outstem', 'with_bias', 'ker_param', 'lambda_set', 'threshold_set', 'metric'};

for i = 1:(numel(varnames)-5)
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

% default of with_bias flag
if(length(varargin)<7 || isempty(varargin{7}))
    param.with_bias = 1;
else
    param.with_bias = varargin{7};
end

% default of ker_param
if(length(varargin)<8 || isempty(varargin{8}) || strcmpi(varargin{8}, 'none'))
    param.ker_param.type = 'corr';
    param.ker_param.scale = NaN;
else
    load(varargin{8})
    param.ker_param = ker_param;
end

% default of lambda_set
if(length(varargin)<9 || isempty(varargin{9}) || strcmpi(varargin{9}, 'none'))
    param.lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 ...
        2.5 3 3.5 4 5 10 15 20];
else
    load(varargin{9})
    param.lambda_set = lambda_set;
end

% default of threshold_set
if(length(varargin)<10 || isempty(varargin{10}) || strcmpi(varargin{10}, 'none'))
    param.threshold_set = [-1:0.1:1];
else
    load(varargin{10})
    param.threshold_set = threshold_set;
end

% default of metric
if(length(varargin)<11 || isempty(varargin{11}))
    param.metric = 'Predictive_COD';
else
    param.metric = varargin{11};
end

end
