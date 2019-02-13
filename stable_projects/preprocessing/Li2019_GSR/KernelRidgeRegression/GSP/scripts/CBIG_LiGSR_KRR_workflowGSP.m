function CBIG_LiGSR_KRR_workflowGSP( data_csv, subject_list, RSFC_file,  y_list, ...
   covariate_list, FD_file, DVARS_file, outdir, outstem, num_test_folds, num_inner_folds, ...
   seed, ker_param_file, lambda_set_file, threshold_set_file )

% CBIG_LiGSR_KRR_workflowGSP( data_csv, subject_list, RSFC_file, y_list, ...
%     covariate_list, outdir, outstem, num_test_folds, num_inner_folds, seed, ...
%     ker_param_file, lambda_set_file, threshold_set_file)
% 
% This function performs the whole kernel ridge regression procedure for
% the Brain Genomics Superstruct Project (GSP) dataset, given a set of
% target traits to be predicted (y_list). It first splits the data into
% cross-validation folds, then read in the traits and covariates, and calls
% "CBIG_KRR_workflow.m" to run kernel ridge regression algorithm. The
% prediction accuracy using the optimal hyperparameters will be saved as
% [outdir '/final_result_' outstem '.mat'].
% 
% Inputs:
%   - data_csv
%     Full path of the CSV file containing behavioral and demographic
%     information from the GSP dataset.
% 
%   - subject_list
%     Full path of the subject ID list. Each line in this list corresponds
%     to one subject.
%  
%   - RSFC_file
%     Full path of the resting-state functional connectivity (RSFC) matrix.
%  
%   - y_list
%     Full path to a text file with all y, i.e. behavioral (or demographic)
%     measures (target measures to be predicted using kernel ridge
%     regression). Each line in this text file corresponds to one
%     behavioral name. The behavioral names should correspond to the
%     headers in "data_csv".
% 
%   - covariate_list
%     A text list of covariate names (e.g. age, sex, FD) that need to be
%     regressed from y (i.e. the measures to be predicted).
%     Each line in this text file corresponds to one covariate name. The
%     covariate names stated in this list, except for 'FD' and 'DVARS',
%     should correspond to the headers in "data_csv".
% 
%   - FD_file (optional)
%     If there is a need to regress FD (framewise displacement) from the
%     behavioral (or demographic) measures, the user should include 'FD' in
%     the "covariate_list". In this case, "FD_file" is the full path of the
%     mean FD of all subjects. The number of lines in "FD_file" should be
%     the same as the number of lines in "subject_list". 
%     If the user does not need to regress FD from y, then the input
%     variable "FD_file" is not required and the user can pass in 'NONE' to
%     the function.
%     If "covariate_list" does not contain FD, this argument will be
%     ignored.
% 
%   - DVARS_file (optional)
%     If there is a need to regress 'DVARS' from the behavioral
%     (demographic) measures, y, the user must include the covariate
%     'DVARS' (or 'DV') in the 'covariate_list'. In this case, "DVARS_file"
%     is the full path of the mean DVARS of all subjects. The number of
%     lines in "DVARS_file" should be the same as the number of lines in
%     "subject_list". 
%     If the user does not need to regress DVARS from y, then the input
%     variable 'DVARS_file' is not required and the user can pass in 'NONE'
%     to the function.
%     If "covariate_list" does not contain DV (or DVARS), this argument
%     will be ignored.
% 
%   - outdir
%     The full path of output directory. A subfolder 
%     [outdir '/randseed_' seed] will be created to save all output files of
%     current random seed.
% 
%   - outstem
%     A string appended to the output file names to specify the output
%     files (y after regression, accuracy files, ...). For example, if
%     outstem = '58behaviors', then the accuracy files will be names as
%     <path_to_file>/acc_23behaviors.mat, 
%     and the final output filename will be 
%     [outdir '/randseed_' seed '/final_result_23behaviors.mat'].
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
%     Full path of the file storing (.mat) the set of threshold used to 
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
fprintf('[GSP workflow]: split families into %d folds.\n', num_test_folds);
CBIG_cross_validation_data_split( subject_list, 'NONE', 'NONE', ...
    'NONE', num_test_folds, seed, fullfile(outdir, ['randseed_' num2str(seed)]), ',' );

%% Read y
% y types
fprintf('[GSP workflow]: read the measures to be predicted.\n')
[y_names, num_y] = CBIG_text2cell(y_list);
for i = 1:num_y
    if(strcmp(y_names{i}, 'Sex'))
        y_types{i} = 'categorical';
    else
        y_types{i} = 'continuous';
    end
end

MenTot_flag = 0;
if(ismember({'avg_MenRot_non0_CORRpc'}, y_names))
    MenTot_flag = 1;
    idx = find(strcmp(y_names,'avg_MenRot_non0_CORRpc')==1);
    y_names = [y_names(1:idx-1) {'MenRot_80_CORRpc' 'MenRot_120_CORRpc' 'MenRot_160_CORRpc'} ...
        y_names(idx+1:end)];
    y_types = [y_types(1:idx-1) {'continuous' 'continuous' 'continuous'} y_types(idx+1:end)];
end

ystem = outstem;
if(~isempty(outstem))
    ystem = ['_' outstem];
end
if(~exist(fullfile(outdir, ['y' ystem '.mat']), 'file'))
    CBIG_read_y_from_csv( {data_csv}, 'Subject_ID', y_names, y_types, ...
        subject_list, fullfile(outdir, ['y' ystem '.mat']), ',' );
    
    if(MenTot_flag == 1)
        y_tmp = load(fullfile(outdir, ['y' ystem '.mat']));
        
        y_tmp.y(:,idx)=mean(y_tmp.y(:,idx:idx+2),2);
        y_tmp.y(:,idx+1:idx+2)=[];
        
        save(fullfile(outdir, ['y' ystem '.mat']), '-struct' ,'y_tmp')
    end
end



%% Read covariates
% covariate types
fprintf('[GSP workflow]: read covariates to be regressed from the measures.\n')
[cov_names, num_cov] = CBIG_text2cell(covariate_list);
for i = 1:num_cov
    if(strcmp(cov_names{i}, 'Sex') || strcmp(cov_names{i}, 'Race_Ethn'))
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
    CBIG_generate_covariates_from_csv( {data_csv}, ...
        'Subject_ID', cov_names, cov_types, subject_list, FD_file, DVARS_file, ...
        fullfile(outdir, ['covariates' cov_stem '.mat']), ',' );
end

%% Call kernel regression workflow utility function
fprintf('[GSP workflow]: call kernel regression workflow ...\n')
sub_fold_file = fullfile(outdir, ['randseed_' num2str(seed)], ...
    ['no_relative_' num2str(num_test_folds) '_fold_sub_list.mat']);
CBIG_KRR_workflow( '', 0, sub_fold_file, fullfile(outdir, ['y' ystem '.mat']), ...
    fullfile(outdir, ['covariates' cov_stem '.mat']), RSFC_file, num_inner_folds, ...
    fullfile(outdir, ['randseed_' num2str(seed)]), outstem, ker_param_file, ...
    lambda_set_file, threshold_set_file);


end

