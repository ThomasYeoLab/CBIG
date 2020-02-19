function CBIG_LiGSR_LRR_workflowHCP( restricted_csv, unrestricted_csv, subject_list, FC_name, ...
    y_name, covariate_list, seed, num_test_folds, num_inner_folds, FD_file, DVARS_file, ...
    gpso_dir, outdir, eval, tree, varargin)

% CBIG_LiGSR_LRR_workflowHCP( restricted_csv, unrestricted_csv, subject_list, FC_name, ...
%     y_name, covariate_list, seed, num_test_folds, num_inner_folds, FD_file, DVARS_file, ...
%     gpso_dir, outdir, eval, tree, varargin)
% 
% This function performs the whole linear ridge regression procedure for
% the Human Connectome Project (GSP) dataset, given a target trait to be
% predicted (y_name). It first splits the data into cross-validation folds,
% then read in the traits and covariates, and calls
% "CBIG_LRR_workflow_1measure.m" to run linear ridge regression algorithm. 
% The prediction accuracy using the optimal hyperparameters will be saved as
% fullfile(outdir, 'results', 'optimal_acc', [y_name '.mat']).
% 
% Inputs:
%   - resitricted_csv
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
%   - FC_name
%     Full path of the resting-state functional connectivity (RSFC) matrix.
% 
%   - y_name
%     A string. The name of the behavior (or trait) to be predicted. The
%     name hould exist as a header in either "restricted_csv" or
%     "unrestricted_csv".
% 
%   - covariate_list
%     Full path to a text file stating all covariate names. Each line
%     corresponds to one covariate name. The covariate names should exist
%     as header in either "restricted_csv" or "unrestricted_csv", except
%     for 'FD' and 'DVARS'.
% 
%   - seed
%     A string or scalar, the random seed used to split the data into
%     training-test cross-validation folds.
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
%   - FD_file (optional) 
%     If there is a need to regress 'FD' from the behavioral (or demographic)
%     measures, y, the user should include 'FD' in the "covariate_list". In
%     this case, "FD_file" is the full path of a text file containing the
%     mean framewise displacement (FD) of all subjects. The number of lines
%     in "FD_file" should be the same as the number of lines in "subject_list". 
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
%     is the full path of a text file containing the mean DVARS of all
%     subjects. The number of lines in "DVARS_file" should be the same as
%     the number of lines in "subject_list". 
%     If the user does not need to regress 'DVARS' from y, then the input
%     variable 'DVARS_file' is not required and the user can pass in 'NONE'
%     to the function.
%     If "covariate_list" does not contain DV (or DVARS), this argument
%     will be ignored.
% 
%   - gpso_dir
%     A string, the full path of the directory to store the cloned git
%     repository for Gaussian process code. Current script needs the
%     Gaussian-Process Surrogate Optimisation (GPSO;
%     https://github.com/jhadida/gpso) package to optimize the objective
%     function. This git repository and its dependency "deck"
%     (https://github.com/jhadida/deck) need to be both cloned into the
%     same folder, which is passed into current script through gpso_dir.
% 
%   - outdir
%     The full path of output directory. A subfolder 
%     [outdir '/randseed_' seed] will be created to save all output files of
%     current random seed.
% 
%   - eval (optional)
%     A string or scalar, the maximal evaluation times of the objective
%     function. Default is 15. If it is too large, the runtime would be too
%     long; if it is too small, the objective function may not be able to
%     reach its optimal value. The users need to consider this trade-off
%     when setting this variable. For more information, please refer to: 
%     https://github.com/jhadida/gpso
% 
%   - tree (optional)
%     A string or scalar, the depth of the partition tree used to explore
%     hyperparameters. Default is 3. For more information, please refer to:
%     https://github.com/jhadida/gpso
% 
%   - varargin (optional)
%     One or more strings. 
%     Some demographic variables (e.g. age) are usually treated as
%     covariates when predicting other behavioral variables (e.g. IQ). In
%     this case, age values are saved together with other covariates in a
%     mat file called 'covariates.mat'.
%     However, sometimes these variables could also be treated as target
%     variables. For example, we want to predict age so that age cannot be
%     included in covariates files anymore. In this case, this script
%     automatically save the covariate file as 'covariate_Age_in_Yrs.mat',
%     where '_Age_in_Yrs' file stem indicates that this covariate file is
%     used to predict age.
%     To achieve this, this script claims a sets of such variable names
%     that could be used as either covariate or target. The default set of
%     such variables include 'Gender', 'FD', 'DVARS', 'Age_in_Yrs', 'Race',
%     'SSAGA_Educ'. If the user have other possible variables, he/she can
%     pass them in through varargin. This script will concatenate varargin
%     with the variables aforementioned.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(ischar(seed))
    seed = str2num(seed);
end
if(ischar(num_test_folds))
    num_test_folds = str2num(num_test_folds);
end
if(ischar(num_inner_folds))
    num_inner_folds = str2num(num_inner_folds);
end

if(~exist('eval', 'var') || isempty(eval))
    eval = 15;
elseif(ischar(eval))
    eval = str2num(eval);
end
if(~exist('tree', 'var') || isempty(tree))
    tree = 3;
elseif(ischar(tree))
    tree = str2num(tree);
end

outstem = y_name;
opt_out = fullfile(outdir, 'results', 'optimal_acc', [outstem '.mat']);

%% Data split
fprintf('[HCP workflow]: split families into %d folds.\n', num_test_folds);
sub_fold_file = fullfile(outdir, ['randseed_' num2str(seed)], ...
    ['no_relative_' num2str(num_test_folds) '_fold_sub_list.mat']);
if(~exist(sub_fold_file, 'file'))
    CBIG_cross_validation_data_split( subject_list, restricted_csv, 'Subject', 'Family_ID', ...
        num_test_folds, seed, fullfile(outdir, ['randseed_' num2str(seed)]), ',' );
end

%% read y
ystem = ['_' y_name];

% y types
fprintf('[HCP workflow]: read the measures to be predicted.\n')
if(strcmp(y_name, 'Gender'))
    y_type = {'categorical'};
else
    y_type = {'continuous'};
end

y_name = {y_name};

if(~exist(fullfile(outdir, 'y', ['y' ystem '.mat']), 'file'))
    CBIG_read_y_from_csv( {restricted_csv, unrestricted_csv}, 'Subject', y_name, y_type, ...
        subject_list, fullfile(outdir, 'y', ['y' ystem '.mat']), ',' );
end

%% read covariates
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
cov_stem = '';
possib_cov = {'Gender', 'FD', 'DVARS', 'Age_in_Yrs', 'Race', 'SSAGA_Educ'};
if(~isempty(varargin))
    possib_cov = [possib_cov varargin];
end
for i = 1:length(possib_cov)
    if(strcmp(outstem, possib_cov{i}))
        cov_stem =['_' possib_cov{i}];
    end
end
if(~exist(fullfile(outdir, 'covariates', ['covariates' cov_stem '.mat']), 'file'))
    CBIG_generate_covariates_from_csv( {restricted_csv, unrestricted_csv}, ...
        'Subject', cov_names, cov_types, subject_list, FD_file, DVARS_file, ...
        fullfile(outdir, 'covariates', ['covariates' cov_stem '.mat']), ',' );
end

%% call general linear ridge regression scripts
if(~exist(opt_out, 'file'))
    load(sub_fold_file)
    params.sub_fold = sub_fold;
    clear sub_fold
    
    load(fullfile(outdir, 'y', ['y' ystem '.mat']))
    params.y = y;
    clear y
    
    load(fullfile(outdir, 'covariates', ['covariates' cov_stem '.mat']))
    params.covariates = covariates;
    clear covariates
    
    load(FC_name)
    params.FC_mat = corr_mat;
    clear corr_mat
    
    params.num_inner_folds = num_inner_folds;
    params.outdir = fullfile(outdir, ['randseed_' num2str(seed)]);
    params.outstem = outstem;
    params.gpso_dir = gpso_dir;
    params.domain = [0.001 0.1; 3 8];
    params.eval = eval;
    params.tree = tree;
    
    [ acc_corr_train, acc_corr_test, mse_test, y_predict] = ...
        CBIG_LRR_workflow_1measure( params );
end

end

