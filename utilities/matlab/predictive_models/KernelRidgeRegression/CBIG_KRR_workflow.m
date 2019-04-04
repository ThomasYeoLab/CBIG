function CBIG_KRR_workflow( setup_file, save_setup, varargin )

% CBIG_KRR_workflow( setup_file, save_setup, varargin )
% 
% This function runs the whole workflow of kernel ridge regression
% algorithm to predict some target variables. All the input filenames and
% hyperparameters can be passed in through the setup_file. Alternatively,
% they can also be passed in one by one through varargin.
% 
% Note: if your target variables consist of both binary and non-binary
% variables, you need to deal with them separately and run this workflow
% twice.
% 
% Inputs:
%   - setup_file
%     Full path of the setup file (.mat storing a structure). If this
%     argument is passed in, all other arguments are not needed because
%     they are assumed to be saved in the setup file. If this argument is
%     not passed in, this function will automatically write out a setup
%     file in the output directory (i.e. varargin{6}) for the convienence
%     of the users to rerun this function. 
%     setup_file has the following fields:
%     1. sub_fold 
%        A structure containing the information of how the data are
%        separated into training and test sets. See the description of
%        varargin{1}.
%     2. y
%        A matrix of the target variables to be predicted. See the
%        description of varargin{2}.
%     3. covariates
%        A matrix of the covariates to be regressed out from y. See the
%        description of varargin{3}.
%     4. feature_mat
%        A matrix of the features used as the independent variable in the
%        prediction. See the description of varargin{4}.
%     5. num_inner_folds
%        A scalar, the number of inner-loop cross-validation folds. See the
%        description of varargin{5}.
%     6. outdir
%        A string of the output directory. See the description of
%        varargin{6}.
%     7. outstem
%        A string appended to the filenames. See the description of
%        varargin{7}.
%     8. ker_param 
%        A structure of all possible kernel parameters. See the description
%        of varargin{8}.
%     9. lambda_set
%        A vector of all possible regularization parameters used for grid
%        search. See the description of varargin{9}.
%     10.threshold_set
%        A vector of all possible thresholds to determine the separation
%        point for binary target variables in the prediction. See the
%        description of varargin{10}.
% 
%   - save_setup
%     A string or a scalar of 0 or 1. If the user passed in 1, then a
%     setup_file will be saved out for the user to rerun if needed. If the
%     user passed in 0, no setup_file will be saved.
%     
%%%%% varargin:
%     In matlab, "varargin" grabs all the parameters passed in, starting
%     from the position of "varargin" till the last input, and store them
%     in cell arrays. In the case of "CBIG_KRR_workflow.m", since
%     "varargin" is the 3rd argument, the 3rd to the last inputs will all
%     be stored in "varargin".
%     Maximum 10 parameters can be passed in through "varargin". If
%     setup_file is not passed in, the user must include the details of the
%     first 7 parameters. The details of each parameter (indicated as a
%     cell array in "varargin") is as stated below:
% 
%   - varargin{1}   (sub_fold_file)
%     Full path of the cross-validation data split file. A structure
%     "sub_fold" is assumed to be saved in this file. 
%     sub_fold(i).fold_index is a #subjects x 1 binary vector. 1 refers to
%     the corresponding subject is a test subject in the i-th test fold. 0
%     refers to the corresponding subject is a training subject in the i-th
%     test fold.
%     If varargin{1} is not passed in, then "sub_fold" structure is assumed
%     to be stored in "setup_file".
% 
%   - varargin{2}   (y_file)
%     Full path of the file storing the original target measures (before
%     regressing out covariates), y, to be used for prediction.
%     A #subjects x #MeasuresToPredict matrix "y" is assumed to be saved in
%     this file.
%     If this file is not passed in, then the matrix "y" is assumed to be
%     stored in "setup_file".
% 
%   - varargin{3}   (covariate_file)
%     Full path of the file storing the covariates that need to be
%     regressed from "y" (age, sex, ...). A #subjects x #regressors matrix
%     "covariates" is assumed to be saved in this file.
%     If this file is not passed in, then the matrix "covariates" is
%     assumed to be stored in "setup_file".
%  
%   - varargin{4}   (feature_file)
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
%   - varargin{5}   (num_inner_folds)
%     A string or a scalar. The number of inner-loop folds within each
%     training fold for hyperparameter selection.
%     If this argument is not passed in, it is assumed to be saved in
%     "setup_file".
% 
%   - varargin{6}   (outdir)
%     Full path of the output directory. If this argument is not passed in,
%     then this string is assumed to be saved in "setup_file".
% 
%   - varargin{7}   (outstem)
%     A string appended to the filename to specify the output files (y
%     after regression, accuracy files, ...). For example, if outstem =
%     '58behaviors', then the accuracy files will be names as
%     <path_to_file>/acc_58behaviors.mat, and the final output filename
%     will be [outdir '/final_result_58behaviors.mat'].
% 
%   - varargin{8}   (ker_param_file)
%     Full path of the kernel parameter file (.mat). A structure "ker_param" 
%     is assumed to be saved in this file.
%     "ker_param" is a K x 1 structure with two fields: type and scale. K
%     denotes the number of kernels.
%     ker_param(k).type is a string that specifies the type of k-th kernel.
%                       Choose from
%                       'corr'        - Pearson's correlation;
%                       'Gaussian'    - Gaussian kernel;
%                       'Exponential' - exponential kernel.
%     ker_param(k).scale is a scalar specifying the scale of k-th kernel
%     (for Gaussian kernel or exponential kernel). If ker_param(k).type == 'corr',
%     ker_param(k).scale = NaN.
%     If this argument is not passed in, and "ker_param" also does NOT
%     exist in "setup_file", then ker_param will be set as default:
%     ker_param.type = 'corr';
%     ker_param.scale = NaN.
% 
%   - varargin{9}   (lambda_set_file)
%     Full path of the regularization parameter file (.mat). A vector 
%     "lambda_set" is assumed to be saved in this file.
%     "lambda_set" is a vector of numbers for grid search of lambda (the
%     regularization parameter). If this file is not passed in and also
%     does NOT exist in "setup_file", or "lambda_set" is 'NONE', it will be
%     set as default:
%     [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
%        5 10 15 20 30 40 50 60 70 80 100 150 200 300 500 700 1000 10000 100000 1000000]
% 
%  - varargin{10}   (threshold_set_file)
%    Full path of the file (.mat) storing the set of threshold used to 
%    binarize the predicted score when the original y is binary. A vector
%    "threshold_set" is assumed to be saved in this file.
%    "threshold_set" is a vector used for grid search of optimal
%    "threshold". If this file is not passed in and also does NOT exist in
%    "setup_file", or "threshold_set" is 'NONE', it will be set as default: 
%    [-1:0.1:1].
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Input arguments
if(~isempty(setup_file))
    param = load(setup_file);
else
    param = CBIG_KRRworkflow_parseInput(varargin);
    if(save_setup==1 || strcmp(save_setup, '1'))
        stem = param.outstem;
        if(~isempty(param.outstem))
            stem = ['_' stem];
        end
        save(fullfile(param.outdir, ['setup' stem '.mat']), '-struct', 'param')
    end
end

num_test_folds = length(param.sub_fold);

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

%% step 1. Regress covariates from y
fprintf('# step 1: regress covariates from y for each fold.\n')
CBIG_crossvalid_regress_covariates_from_y( ...
    param.y, param.covariates, param.sub_fold, param.outdir, param.outstem);

%% step 2. Generate kernels
fprintf('# Step 2: generate kernels.\n')
CBIG_KRR_generate_kernels( param.feature_mat, param.sub_fold, param.outdir, param.ker_param )

%% step 3. Inner-loop cross-validation
fprintf('# Step 3: inner-loop cross-validation.\n')
for test_fold = 1:num_test_folds
    CBIG_KRR_innerloop_cv_allparams( test_fold, param.sub_fold, param.num_inner_folds, ...
        param.outdir, param.outstem, param.ker_param, param.lambda_set, param.threshold_set );
end

%% step 4. Test cross-validation
fprintf('# Step 4: training-test cross-validation.\n')
for test_fold = 1:num_test_folds
    CBIG_KRR_test_cv_allparams( test_fold, param.sub_fold, ...
        param.outdir, param.outstem, param.ker_param, param.lambda_set, param.threshold_set );
end

%% step 5. Select optima
fprintf('# Step 5: select optimal hyperparameters and obtain the test accuracy.\n')
CBIG_KRR_pick_optima( param.sub_fold, param.outdir, param.outstem, param.bin_flag, ...
    param.ker_param, param.lambda_set, param.threshold_set );






function param = CBIG_KRRworkflow_parseInput(var)

varnames = {'sub_fold', 'y', 'covariates', 'feature_mat', 'num_inner_folds', 'outdir', ...
    'outstem', 'ker_param', 'lambda_set', 'threshold_set'};

for i = 1:(numel(varnames)-3)
    if(length(var)<i || isempty(var{i}))
        error('Variable ''%s'' is needed.', varnames{i})
    else
        if(i~=5 && i~=6 && i~=7)
            % sub_fold, y, covariates, feature_mat are loaded from .mat files.
            currvar = load(var{i});
            names = fieldnames(currvar);
            if(any(strcmp(names, varnames{i})))
                currname = varnames{i};
            else
                currname = names{1};
            end
            param.(varnames{i}) = currvar.(currname);
        else
            % num_inner_folds or outdir are directly set as fields of param
            param.(varnames{i}) = var{i};
        end
    end
end

% default of ker_param
if(length(var)<8 || isempty(var{8}) || strcmpi(var{8}, 'none'))
    param.ker_param.type = 'corr';
    param.ker_param.scale = NaN;
else
    load(var{8})
    param.ker_param = ker_param;
end

% default of lambda_set
if(length(var)<9 || isempty(var{9}) || strcmpi(var{9}, 'none'))
    param.lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 ...
        2.5 3 3.5 4 5 10 15 20 30 40 50 60 70 80 100 150 200 300 500 700 1000 10000 100000 1000000];
else
    load(var{9})
    param.lambda_set = lambda_set;
end

% default of threshold_set
if(length(var)<10 || isempty(var{10}) || strcmpi(var{10}, 'none'))
    param.threshold_set = [-1:0.1:1];
else
    load(var{10})
    param.threshold_set = threshold_set;
end



