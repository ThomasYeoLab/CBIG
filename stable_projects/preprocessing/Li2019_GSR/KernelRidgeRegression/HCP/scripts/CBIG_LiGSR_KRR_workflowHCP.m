function CBIG_LiGSR_KRR_workflowHCP( restricted_csv, unrestricted_csv, subject_list, RSFC_file, ...
    y_list, covariate_list, FD_file, DVARS_file, outdir, outstem, num_test_folds, num_inner_folds, ...
    seed, ker_param_file, lambda_set_file, threshold_set_file )

% CBIG_LiGSR_KRR_workflowHCP( restricted_csv, unrestricted_csv, subject_list, RSFC_file, ...
%     y_list, covariate_list, outdir, outstem, num_test_folds, num_inner_folds, seed, ...
%     ker_param_file, lambda_set_file, threshold_set_file)
% 
% This function performs the whole kernel ridge regression procedure for
% the Human Connectome Project (HCP) dataset, given a set of target traits
% to be predicted (y_list). It first splits the data into cross-validation
% folds, then read in the traits and covariates, and calls
% "CBIG_KRR_workflow.m" to run kernel ridge regression algorithm. The
% prediction accuracy using the optimal hyperparameters will be saved as
% [outdir '/final_result_' outstem '.mat'].
% 
% Inputs:
%   - restricted_csv
%     Full path of the restricted CSV file downloaded from the HCP website.
% 
%   - unrestricted_csv
%     Full path of the unrestricted CSV file downloaded from the HCP
%     website.
% 
%   - subject_list
%     Full path of the subject ID list. Each line in this list corresponds
%     to one subject ID.
%  
%   - RSFC_file
%     Full path of the resting-state functional connectivity (RSFC) matrix.
%  
%   - y_list
%     Full path to a text file with all behavioral (or demographic)
%     measures (measures to be predicted using kernel ridge regression). 
%     Each line corresponds to one behavioral name. The behavioral names
%     should exist as a header in either "restricted_csv"
%     or "unrestricted_csv".
% 
%   - covariate_list
%     Full path to a text file stating all covariate names. Each line
%     corresponds to one covariate name. The covariate names should exist
%     as header in either "restricted_csv" or "unrestricted_csv", except
%     for 'FD' and 'DVARS'.
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
%   - num_test_folds
%     A string or scalar, the number of training-test cross-validation
%     folds.
% 
%   - num_inner_folds
%     A string or scalar.
%     To select optimal hyperparameters, each training fold will be split
%     randomly into "num_inner_folds" inner-loop cross-validation folds.
% 
%   - seed
%     A string or scalar, the random seed used to split the data into
%     training-test cross-validation folds.
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
%        5 10 15 20 30 40 50 60 70 80 100 150 200 300 500 700 1000 10000 100000 1000000]
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
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% setting up
if(ischar(num_test_folds))
    num_test_folds = str2double(num_test_folds);
end

if(ischar(num_inner_folds))
    num_inner_folds = str2double(num_inner_folds);
end

if(ischar(seed))
    seed = str2double(seed);
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
fprintf('[HCP workflow]: split families into %d folds.\n', num_test_folds);
CBIG_cross_validation_data_split( subject_list, restricted_csv, 'Subject', ...
    'Family_ID', num_test_folds, seed, fullfile(outdir, ['randseed_' num2str(seed)]), ',' );

%% Read y
% y types
fprintf('[HCP workflow]: read the measures to be predicted.\n')
[y_names, num_y] = CBIG_text2cell(y_list);
for i = 1:num_y
    if(strcmp(y_names{i}, 'Gender'))
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
    CBIG_read_y_from_csv( {restricted_csv, unrestricted_csv}, 'Subject', y_names, y_types, ...
        subject_list, fullfile(outdir, ['y' ystem '.mat']), ',' );
end

%% Read covariates
% covariate types
fprintf('[HCP workflow]: read covariates to be regressed from the measures.\n')
[cov_names, num_cov] = CBIG_text2cell(covariate_list);
for i = 1:num_cov
    if(strcmp(cov_names{i}, 'Gender') || strcmp(cov_names{i}, 'Race'))
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
    CBIG_generate_covariates_from_csv( {restricted_csv, unrestricted_csv}, ...
        'Subject', cov_names, cov_types, subject_list, FD_file, DVARS_file, ...
        fullfile(outdir, ['covariates' cov_stem '.mat']), ',' );
end

%% Call kernel regression workflow utility function
fprintf('[HCP workflow]: call kernel regression workflow ...\n')
sub_fold_file = fullfile(outdir, ['randseed_' num2str(seed)], ...
    ['no_relative_' num2str(num_test_folds) '_fold_sub_list.mat']);
CBIG_KRR_workflow( '', 0, sub_fold_file, fullfile(outdir, ['y' ystem '.mat']), ...
    fullfile(outdir, ['covariates' cov_stem '.mat']), RSFC_file, num_inner_folds, ...
    fullfile(outdir, ['randseed_' num2str(seed)]), outstem, ker_param_file, ...
    lambda_set_file, threshold_set_file);


end

