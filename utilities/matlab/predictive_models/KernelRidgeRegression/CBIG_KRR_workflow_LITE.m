function CBIG_KRR_workflow_LITE( setup_param, save_setup, sub_fold_file, y_file, ...
    covariate_file, feature_file, num_inner_folds, outdir, outstem, varargin )

% CBIG_KRR_workflow_LITE( setup_param, save_setup, sub_fold_file, y_file, ...
%    covariate_file, feature_file, num_inner_folds, outdir, outstem, varargin )
% 
% This function runs the whole workflow of kernel ridge regression
% algorithm to predict some target variables. This is a lite version of
% `CBIG_KRR_workflow.m`. In this version, the kernels are not computed and
% saved separately for each inner-loop fold and each outer training-test
% fold. Instead, the kernel is only computed once across all the subjects,
% and for each fold, the scripts will grab the corresponding entries on the
% fly.
% 
% All the input filenames and hyperparameters can be passed in through the
% setup_param. Alternatively,
% they can also be passed in one by one through the variables after
% 'save_setup'. Optional input variables can be passed in through varargin
% by specifying the 'name' and 'label'. An example would be as follows:
%
% CBIG_KRR_workflow_LITE( '', 1, sub_fold_file, y_file, ...
%    covariate_file, feature_file, num_inner_folds, outdir, outstem, 'with_bias',...
%     0, 'lambda_set_file', 'xxx/xxx/lambda_set_file.mat')
% 
% Note: if your target variables consist of both binary and non-binary
% variables, you need to deal with them separately and run this workflow
% twice.
% 
% Inputs:
%   - setup_param
%     Setup parameter (a structure) or full path of the setup parameter
%     (.mat storing a structure). If this argument is passed in, all other
%     arguments are not needed because they are assumed to be saved in the
%     setup param. If this argument is not passed in, this function will
%     automatically write out a setup file in the output directory
%     (i.e. varargin{6}) for the convienence of the users to rerun this function.
%     setup_param has the following fields:
%     1. sub_fold 
%        A structure containing the information of how the data are
%        separated into training and test sets. See the description of
%        `sub_fold_file` in the Compulsory and Optional Variables section.
%     2. y
%        A matrix of the target variables to be predicted. See the
%        description of `y_file` in the Compulsory and Optional Variables section
%     3. covariates
%        A matrix of the covariates to be regressed out from y. See the
%        description of `covariate_file` in the Compulsory and Optional Variables section.
%     4. feature_mat
%        A matrix of the features used as the independent variable in the
%        prediction. See the description of See the description of `feature_file` 
%        in the Compulsory and Optional Variables section.
%     5. num_inner_folds
%        A scalar, the number of inner-loop cross-validation folds. See the
%        description of `num_inner_folds` in the Compulsory and Optional
%        Variables section.
%     6. outdir
%        A string of the output directory. See the description of
%        `outdir` in the Compulsory and Optional Variables section.
%     7. outstem
%        A string appended to the filenames. See the description of
%        `outstem` in the Compulsory and Optional Variables section.
%     8. with_bias
%        A string (choose from '0' or '1') or a scalar (choose from 0 or 1).
%        See the description of `with_bias` in the Compulsory and Optional 
%        Variables section
%     9. ker_param 
%        A structure of all possible kernel parameters. See the description
%        of `ker_param_file` in the Compulsory and Optional Variables section
%     10.lambda_set
%        A vector of all possible regularization parameters used for grid
%        search. See the description of `lambda_set_file` in the Compulsory
%        and Optional Variables section.
%     11.threshold_set
%        A vector of all possible thresholds to determine the separation
%        point for binary target variables in the prediction. See the
%        description of `threshold_set_file` in the Compulsory and 
%        Optional Variables section.
%     12.metric
%        A string indicating the metric used to define prediction loss. See the
%        description of `metric` in the Compulsory and Optional Variables section.
% 
%   - save_setup
%     A string or a scalar of 0 or 1. If the user passed in 1, then a
%     setup_param will be saved out for the user to rerun if needed. If the
%     user passed in 0, no setup_param will be saved.
%     
%%%%% Compulsory and Optional Variables if `setup_param` is not passed in:
%
%     If setup_param is not passed in, the user must include the details of the
%     first 7 parameters (ie from `sub_fold_file` to `outstem). The details of
%     each parameter is as stated below:
%
%     In matlab, "varargin" grabs all the parameters passed in, starting
%     from the position of "varargin" till the last input, and store them
%     in cell arrays. In the case of "CBIG_KRR_workflow.m", since
%     "varargin" is the 10th argument, the 10th to the last inputs will all
%     be stored in "varargin". This "varargin" will handle optional input
%     variables. All the user needs to do is, for example if they want to
%     specify 'with_bias' to be 0, they just need to pass in 
%     (...,outstem,'with_bias',0,...), ie specify the name of the variable
%     first and then what you want the value of the variable to be.
%
%     Compulsory Variables
%
%   - sub_fold_file
%     Full path of the cross-validation data split file. A structure
%     "sub_fold" is assumed to be saved in this file. 
%     sub_fold(i).fold_index is a #subjects x 1 binary vector. 1 refers to
%     the corresponding subject is a test subject in the i-th test fold. 0
%     refers to the corresponding subject is a training subject in the i-th
%     test fold.
%     `sub_fold_file` is not passed in, then "sub_fold" structure is assumed
%     to be stored in "setup_param".
% 
%   - y_file
%     Full path of the file storing the original target measures (before
%     regressing out covariates), y, to be used for prediction.
%     A #subjects x #MeasuresToPredict matrix "y" is assumed to be saved in
%     this file.
%     If this file is not passed in, then the matrix "y" is assumed to be
%     stored in "setup_param".
% 
%   - covariate_file
%     Full path of the file storing the covariates that need to be
%     regressed from "y" (age, sex, ...). A #subjects x #regressors matrix
%     "covariates" is assumed to be saved in this file.
%     If this file is not passed in, then the matrix "covariates" is
%     assumed to be stored in "setup_param".
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
%     assumed to be stored in "setup_param".
% 
%   - num_inner_folds
%     A string or a scalar. The number of inner-loop folds within each
%     training fold for hyperparameter selection.
%     If this argument is not passed in, it is assumed to be saved in
%     "setup_param".
% 
%   - outdir
%     Full path of the output directory. If this argument is not passed in,
%     then this string is assumed to be saved in "setup_param".
% 
%   - outstem
%     A string appended to the filename to specify the output files (y
%     after regression, accuracy files, ...). For example, if outstem =
%     '58behaviors', then the accuracy files will be names as
%     <path_to_file>/acc_58behaviors.mat, and the final output filename
%     will be [outdir '/final_result_58behaviors.mat'].
%
%     Optional Variables
%   
%   - (...,'with_bias',with_bias,...)
%     A string (choose from '0' or '1') or a scalar (choose from 0 or 1).
%     - with_bias = 0 (or '0') means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 (or '1') means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
%     If not passed in, the default value is 1, i.e. there will be a bias
%     term to be estimated.
%
%   - (...,'save_kernel',save_kernel,...)
%     A string (choose from '0' or '1') or a scalar (choose from 0 or 1).
%     - save_kernel = 0 (or '0') means the algorithm is do not save kernel
%     into files.
%     - save_kernel = 1 (or '1') means the algorithm is save kernel into
%     files.
%     If not passed in, the default value is 1, i.e. kernel will be saved
%     into files.
% 
%   - (...,'ker_param_file',ker_param_file,...)
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
%     exist in "setup_param", then ker_param will be set as default:
%     ker_param.type = 'corr';
%     ker_param.scale = NaN.
% 
%   - (...,'lambda_set_file',lambda_set_file, ...)
%     Full path of the regularization parameter file (.mat). A vector 
%     "lambda_set" is assumed to be saved in this file.
%     "lambda_set" is a vector of numbers for grid search of lambda (the
%     regularization parameter). If this file is not passed in and also
%     does NOT exist in "setup_param", it will be
%     set as default:
%     [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
%        5 10 15 20]
% 
%   - (...,'threshold_set_file',threshold_set_file,...)
%     Full path of the file (.mat) storing the set of threshold used to 
%     binarize the predicted score when the original y is binary. A vector
%     "threshold_set" is assumed to be saved in this file.
%     "threshold_set" is a vector used for grid search of optimal
%     "threshold". If this file is not passed in and also does NOT exist in
%     "setup_param", it will be set as default: 
%     [-1:0.1:1].
%
%   - (...,'metric',metric,...)
%     A string indicating the metric used to define prediction loss. The
%     loss is used to choose hyperparameters.
%     Choose from:
%       'corr'              - Pearson's correlation;
%       'COD'               - Coefficient of determination (R squared)
%       'predictive_COD'    - Predictive coefficient of determination. ---
%                             to do--------------
%       'MAE'               - mean absolute error
%       'MAE_norm'          - mean absolute error divided by the standard
%                             derivation of the target variable of the traning set
%       'MSE'               - mean squared error
%       'MSE_norm'          - mean squared error divided by the variance
%                             of the target variable of the traning set
% 
%  Outputs:
%   Three subfolds and one file will be generated in the output directory: 
%   'y', 'innerloop_cv', 'test_cv', 'final_result_<outstem>.mat'.
%
%    In 'y', this folder contains the original y values and the y values 
%            after the covariates have been regressed for each fold.
% 
%    In 'innerloop_cv', this folder contains the accuracy, loss and predicted 
%            y value for each lambda value for the innerloop. The results 
%            are saved in a separate mat file for each fold.
% 
%    In 'test_cv', this folder contains the accuracy, loss and predicted y 
%            value for each lambda value for the outerloop. The results are 
%            saved in a separate mat file for each fold.
%
%   In 'final_result_<outstem>.mat', it contains final accuracies for the test 
%            folds. A mat file with the following fields are generated:   
%   - optimal_acc 
%     Test accuracy for each fold (given in correlation).
% 
%   - optimal_kernel 
%     A #outerfolds x #behaviors struct containing the kernel type selected 
%     for each test fold.
% 
%   - optimal_lambda 
%     A #outerfolds x #behaviors matrix containing lambda selected for each 
%     test fold.
% 
%   - optimal_threshold 
%     A #outerfolds x #behaviors matrix containing the threshold selected 
%     for each test fold.
% 
%   - y_predict_concat 
%     A #subjects x #behaviors matrix of predicted target values.
% 
%   - optimal_stats 
%     A cell array storing the accuracies of each possible accuracy metric 
%     (eg. corr, MAE, etc). Each cell array is #outerfolds x #behaviors.
%
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Input arguments
if(~isempty(setup_param))
    if isStringScalar(setup_param) || ischar(setup_param)
        param = load(setup_param);
    elseif isstruct(setup_param)
        param = setup_param;
    else
        error('setup_param is not a string or struct')
    end
else
     param_field_names = {'sub_fold', 'y','covariates','feature_mat',...
        'num_inner_folds','outdir','outstem', 'with_bias', 'ker_param',...
        'lambda_set','threshold_set', 'metric', 'save_kernel'};
    compulsory_variable = {'sub_fold_file', 'y_file','covariate_file','feature_file',...
        'num_inner_folds','outdir','outstem'};
    for check_var = 1:length(compulsory_variable)
        bool_var = exist(compulsory_variable{check_var},'var');
        if bool_var == 0
            error([compulsory_variable{check_var} ' is compulsory but not passed in'])
        else
            if check_var < 5
                currvar = load(eval(compulsory_variable{check_var}));
                names = fieldnames(currvar);
                currname = names{1};
                param.(param_field_names{check_var}) = currvar.(currname);
            else
                param.(param_field_names{check_var}) = eval(compulsory_variable{check_var});
            end
        end
    end
    pnames = { 'with_bias' 'ker_param_file' 'lambda_set_file' 'threshold_set_file'...
        'metric' 'save_kernel'};
    dflts =  {1 [] [] [] 'predictive_COD' '1'};
    
    [with_bias, ker_param_file, lambda_set_file, threshold_set_file, metric, save_kernel] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    param.(param_field_names{8}) = with_bias;
    
    if isempty(ker_param_file)
        ker_param.type = 'corr';
        ker_param.scale = NaN;
        param.(param_field_names{9}) = ker_param;
    else
        currvar = load(ker_param_file);
        names = fieldnames(currvar);
        currname = names{1};
        param.(param_field_names{9}) = currvar.(currname);
    end
    
    if isempty(lambda_set_file)
        param.(param_field_names{10}) = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 ...
        2.5 3 3.5 4 5 10 15 20];
    else
        currvar = load(lambda_set_file);
        names = fieldnames(currvar);
        currname = names{1};
        param.(param_field_names{10}) = currvar.(currname);
    end
    
    if isempty(threshold_set_file)
        param.(param_field_names{11}) = [-1:0.1:1];
    else
        currvar = load(threshold_set_file);
        names = fieldnames(currvar);
        currname = names{1};
        param.(param_field_names{11}) = currvar.(currname);
    end
    
    param.(param_field_names{12}) = metric;

    if isempty(threshold_set_file)
        param.(param_field_names{13}) = 1;
    else
        param.(param_field_names{13}) = save_kernel;
    end

    if(save_setup==1 || strcmp(save_setup, '1'))
        stem = param.outstem;
        if(~isempty(param.outstem))
            stem = ['_' stem];
        end
        save(fullfile(param.outdir, ['setup' stem '.mat']), '-struct', 'param')
    end
end

if(ischar(param.with_bias))
    param.with_bias = str2num(param.with_bias);
end

if(isfield(param, 'save_kernel'))
    if(ischar(param.save_kernel))
        param.save_kernel = str2num(param.save_kernel);
    end
else
    param.save_kernel = 1;
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
kernel = CBIG_KRR_generate_kernels_LITE( param.feature_mat, param.outdir, param.save_kernel, param.ker_param )

%% step 3. Inner-loop cross-validation
fprintf('# Step 3: inner-loop cross-validation.\n')
for test_fold = 1:num_test_folds
    CBIG_KRR_innerloop_cv_allparams_LITE( test_fold, param.sub_fold, param.num_inner_folds, ...
        param.outdir, param.outstem, param.with_bias, kernel, param.ker_param, param.lambda_set, ...
        param.threshold_set, param.metric );
end

%% step 4. Test cross-validation
fprintf('# Step 4: training-test cross-validation.\n')
for test_fold = 1:num_test_folds
    CBIG_KRR_test_cv_allparams_LITE( test_fold, param.sub_fold, param.outdir, param.outstem, ...
        param.with_bias, kernel, param.ker_param, param.lambda_set, param.threshold_set );
end

%% step 5. Select optima
fprintf('# Step 5: select optimal hyperparameters and obtain the test accuracy.\n')
CBIG_KRR_pick_optima( param.sub_fold, param.outdir, param.outstem, param.bin_flag, ...
    param.ker_param, param.lambda_set, param.threshold_set );


