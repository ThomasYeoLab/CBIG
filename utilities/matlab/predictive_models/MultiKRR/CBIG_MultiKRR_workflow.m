function CBIG_MultiKRR_workflow( setup_file, save_setup, sub_fold_file, ...
    y_file, covariate_file, feature_mat_file, num_inner_folds, outdir,...
    outstem, varargin)

% CBIG_MultiKRR_workflow( setup_file, save_setup, sub_fold_file, ...
%       y_file, covariate_file, feature_file, num_inner_folds, outdir,...
%       outstem, varargin)
% 
% This function runs the whole workflow of Multi kernel ridge regression
% algorithm to predict some target variables.The work flow includes,
% regressing out covariates, generating the individual kernels and finally
% consolidating the kernels and passing it through the multi KRR model. All
% the input filenames and hyperparameters can be passed in through the
% setup_file. Alternatively,they can also be passed in one by one through
% the compulsory input variables from `sub_fold_file` to `outstem` and also
% optional input variables into `varargin`. To pass in the optional input
% variables, all tehe user needs to do is to pass in the name of the
% optional variable (eg. 'with_bias') and then its corresponding value (eg.
% 0). More concretely, the user can do this CBIG_MultiKRR_workflow('', 1,
% sub_fold_file,...., outstem, 'with_bias', 0).
% 
% Note: if your target variables consist of both binary and non-binary
% variables, you need to deal with them separately and run this workflow
% twice.
% 
% Inputs:
%   - setup_file
%     Full path of the setup file (.mat storing the following variables). If this
%     argument is passed in, varargins are not needed because
%     they are assumed to be saved in the setup file. If this argument is
%     not passed in, this function will consolidate varargins and write out
%     a setup file in the output directory (i.e. varargin{6}) if the user
%     chooses to do so (by setting save_setup to 1), for the convienence of
%     the users to rerun this function.
%     In summary, setup_file has the following variables:
%     1. sub_fold 
%        A structure containing the information of how the data are
%        separated into training and test sets. See the description of
%        `sub_fold_file` in `Compulsory and Optional Variables section`.
%     2. y
%        A matrix of the target variables to be predicted. See the
%        description of `y_file` in `Compulsory and Optional Variables section`.
%     3. covariates
%        A matrix of the covariates to be regressed out from y. See the
%        description of `covariate_file` in `Compulsory and Optional Variables section`.
%     4. feature_mat
%        A cell array of feature matrices used as the independent variable
%        in the  prediction. See the description of `feature_mat_file` in
%        `Compulsory and Optional Variables section`.
%     5. num_inner_folds
%        A scalar, the number of inner-loop cross-validation folds. See the
%        description of `num_inner_folds` in `Compulsory and Optional Variables section`.
%     6. outdir
%        A string of the output directory. See the description of
%        `outdir` in `Compulsory and Optional Variables section`.
%     7. outstem
%        A string appended to the filenames. See the description of
%        `outstem` in `Compulsory and Optional Variables section`.
%     8. with_bias
%        A string (choose from '0' or '1') or a scalar (choose from 0 or 1).
%        See the description of `with_bias` in `Compulsory and Optional Variables section`.
%     9. ker_param 
%        A structure of all possible kernel parameters. See the description
%        of `ker_param` in `Compulsory and Optional Variables section`.
%     10.threshold
%        A vector of all possible thresholds to determine the separation
%        point for binary target variables in the prediction. See the
%        description of `threshold` in `Compulsory and Optional Variables section`.
%     11.group_kernel
%        A cell array of how the user would like to group the various
%        kernels for multi KRR. See the description of `group_kernel_file`
%        in `Compulsory and Optional Variables section`
%     12.domain
%        A 1 x 2 vector specifying the search domain for the
%        hyperparameters of the multi KRR model.The first element of the  
%        vector stores the lower bound and the second element stores the
%        upper bound of the search domain. See the description in `domain_file`
%        in `Compulsory and Optional Variables section`.
%     13. acc_metric
%        A string stating which accuracy metric to be used. See the description in
%        `acc_metric` in the `Compulsory and Optional Variables section`.
% 
%   - save_setup
%     A string or a scalar of 0 or 1. If the user passed in 1, then a
%     setup_file will be saved out for the user to rerun if needed. If the
%     user passed in 0, no setup_file will be saved.
%     
%%%%% Compulsory and Optional Variables if`setup_file` is not passed in:
%
%     If setup_file is not passed in, the user must include the details of the
%     first 7 parameters (ie from `sub_fold_file` to `outstem. The details of
%     each parameter is as stated below:
%
%     "varargin" will handle option variables. "varargin"  grabs all the
%     parameters passed in, starting from the position of "varargin"
%     till the last input, and store them
%     in cell arrays. In the case of "CBIG_MultiKRR_workflow.m", since
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
%     If sub_fold_file is not passed in, then "sub_fold" structure is assumed
%     to be stored in "setup_file".
% 
%   - y_file
%     Full path of the file storing the original target measures (before
%     regressing out covariates), y, to be used for prediction.
%     A #subjects x #MeasuresToPredict matrix "y" is assumed to be saved in
%     this file.
%     If this file is not passed in, then the matrix "y" is assumed to be
%     stored in "setup_file".
% 
%   - covariate_file
%     Full path of the file storing the covariates that need to be
%     regressed from "y" (age, sex, ...). A #subjects x #regressors matrix
%     "covariates" is assumed to be saved in this file.
%     If this file is not passed in, then the matrix "covariates" is
%     assumed to be stored in "setup_file".
%  
%   - feature_mat_file
%     Full path of the file storing a a cell array. The cell array contains
%     the full path of all the feature files that are required to
%     calculate the kernels. For each feature file, a matrix "feature_mat" 
%     is assumed to be saved in this file. "feature_mat" can be a 2-D matrix
%      with dimension of #features x #subjects or a 3-D
%     matrix with dimension of #ROIs1 x #ROIs2 x #subjects. If
%     "feature_mat" is 3-D, it is a connectivity matrix between two sets of
%     ROIs, and only the lower-triangular off-diagonal entries will be
%     considered as features because the connectivity matrix is symmetric.
%     If this file is not passed in, then the matrix "feature_mat" is
%     assumed to be stored in "setup_file".
% 
%   - num_inner_folds
%     A string or a scalar. The number of inner-loop folds within each
%     training fold for hyperparameter selection.
%     If this argument is not passed in, it is assumed to be saved in
%     "setup_file".
% 
%   - outdir
%     Full path of the output directory. If this argument is not passed in,
%     then this string is assumed to be saved in "setup_file".
% 
%   - outstem
%     A string appended to the filename to specify the output files (y
%     after regression, accuracy files, ...). For example, if outstem =
%     '58behaviors', then the accuracy files will be named as
%     <path_to_file>/acc_58behaviors.mat, and the final output filename
%     will be [outdir '/final_result_58behaviors.mat'].
%
%   Optional Variables
%
%   - (...'with_bias', with_bias,...) 
%     A string (choose from '0' or '1') or a scalar (choose from 0 or 1).
%     - with_bias = 0 (or '0') means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 (or '1') means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
%     If not passed in, the default value is 1, i.e. there will be a bias
%     term to be estimated.
% 
%   - (...'ker_param_file', ker_param_file,...) 
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
%   - (...,'threshold',threshold,...) 
%     A string or scalar specifying the threshold used to 
%     binarize the predicted score when the original y is binary. A scalar
%     "threshold" is assumed to be saved in this file. If this file is not
%     passed in and also does NOT exist in
%     "setup_file", or "threshold" is 'NONE', it will be set as default: 
%     0.5.
%
%   - (...,'group_kernel_file',group_kernel_file,...) 
%     Full path of the file (.mat) storing how the user would like to group
%     the various kernels for multi KRR. The file would store a
%     #kernel_groups cell array. Each element in the cell array contains   
%     a vector of size 1 x #kernels_in_the_group. For example if 4 kernels
%     are to be used to train the multi kernel model, but they can be
%     grouped such that 3 of the kernels uses the same hyperparameter
%     lambda for regularization, then I can pass in the 'group_kernel' as
%     {[1,3,4],[2]}. This means that kernels 1, 3 and 4 are regularized
%     using the same lambda while kernel 2 is regularized using a different
%     lambda. Do note that kernel 1 would correspond to the first feature
%     mat in 'feature_file (varargin{4}) so on and so forth.
%     If this file is not passed in and also does NOT exist in "setup_file"
%     or group_kernel is empty ,i.e.
%     '', then it is assumed that each kernel will be regularised
%     independently.  
%
%   - (...,'domain_file',domain_file,...) 
%     Full path of the file (.mat) storing the domain used to 
%     specify the search boundary of the hyperparameter. A 1 x 2 vector
%     "domain" is assumed to be saved in this file.If this file is not
%     passed in and also does NOT exist in
%     "setup_file", or "domain" is 'NONE', it will be set as default: 
%     [0 20]. The first element of the domain stores the lower bound and
%     the second element stores the upper bound
%
%   - (...,'acc_metric',acc_metric,...) 
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
%    Outputs:
%    The optimal test prediction accuracies ('acc') would be stored in the folder
%    [outdir '/final_result_' stem '.mat']. acc would be a cell array of
%    size K where K is the number of kernel types (see ker_param
%    description). Each cell array would contain a #test_fold by
%    #target_behaviors matrix of optimal test prediction accuracies.
%
%    optimal_stats
%     A #k sized cell array. Each element in the cell array stores a data
%     structure. There are 2 fields in the data structure, 'value' and
%     'description'. 'description' describes which accuracy metric the
%     'value' belongs to. 'value' will be a # folds x #behavioral_measures vector
%     storing the values of the optimal test prediction accuracies
%     corresponding to the specified acuracy metric. Within each field,
%     there will be a total number of 7 entries, each entry corresponding
%     to the possible acc_metrics.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Yan Rui Tan and Ru(by) Kong

%% Input arguments

if(~isempty(setup_file))
    clearvars -except setup_file
    param = load(setup_file);
else
    param_field_names = {'sub_fold', 'y','covariates','feature_mat',...
        'num_inner_folds','outdir','outstem', 'with_bias', 'ker_param',...
        'threshold','group_kernel', 'domain','acc_metric'};
     compulsory_variable = {'sub_fold_file', 'y_file','covariate_file','feature_mat_file',...
        'num_inner_folds','outdir','outstem'};
    for check_var = 1:length(compulsory_variable)
        bool_var = exist(compulsory_variable{check_var},'var');
        if bool_var == 0
            error([compulsory_variable{check_var} ' is compulsory but not passed in'])
        else
            if check_var < 4
                currvar = load(eval(compulsory_variable{check_var}));
                names = fieldnames(currvar);
                currname = names{1};
                param.(param_field_names{check_var}) = currvar.(currname);
            elseif check_var == 4
                currvar = load(eval(compulsory_variable{check_var}));
                names = fieldnames(currvar);
                currname = names{1};
                feature_mat_path = currvar.(currname);
                for ii = 1:length(feature_mat_path)
                    currvar = load(feature_mat_path{ii});
                    names = fieldnames(currvar);
                    currname = names{1};
                    feature_mat{ii} = currvar.(currname);
                end
                param.(param_field_names{check_var}) = feature_mat;
            else
                param.(param_field_names{check_var}) = eval(compulsory_variable{check_var});
            end
        end
    end
    
    pnames = { 'with_bias'  'ker_param_file' 'threshold' 'group_kernel_file'...
        'domain_file' 'acc_metric'};
    dflts =  {1 [] 0.5 [] [] 'predictive_COD'};

    [with_bias, ker_param_file, threshold, group_kernel_file, domain_file,...
        acc_metric] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
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
    
    param.(param_field_names{10}) = threshold;
    
    if isempty(group_kernel_file)
        for i = 1:length(feature_mat_path)
            group_kernel{i} = i;
        end
        param.(param_field_names{11}) = group_kernel;
    else
        currvar = load(group_kernel_file);
        names = fieldnames(currvar);
        currname = names{1};
        param.(param_field_names{11}) = currvar.(currname);
    end
    
    if isempty(domain_file)
        param.(param_field_names{12}) = [0 20];
    else
        currvar = load(domain_file);
        names = fieldnames(currvar);
        currname = names{1};
        param.(param_field_names{12}) = currvar.(currname);
    end
    
    param.(param_field_names{13}) = acc_metric;

    if(save_setup==1 || strcmp(save_setup, '1'))
        stem = param.outstem;
        if(~isempty(param.outstem))
            stem = ['_' stem];
        end
        save(fullfile(param.outdir, ['setup' stem '.mat']), '-struct', 'param')
    end
end

%% Check dimensions of inputs

if(ischar(param.with_bias))
    param.with_bias = str2num(param.with_bias);
end

if(ischar(param.threshold))
    param.threshold = str2num(param.threshold);
end

if(ischar(param.num_inner_folds))
    param.num_inner_folds = str2num(param.num_inner_folds);
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
    param.threshold = NaN;
end

%% step 1. Regress covariates from y
fprintf('# step 1: regress covariates from y for each fold.\n')
CBIG_crossvalid_regress_covariates_from_y( ...
    param.y, param.covariates, param.sub_fold, param.outdir, param.outstem);

%% step 2. Generate kernels
fprintf('# Step 2: generate kernels.\n')
kernel_folders = cell(1,length(param.feature_mat));
for i = 1:length(param.feature_mat)
    feat_mat = param.feature_mat{i};
    kernel_folders{i} = ['Kernel_' num2str(i)];
    kernel_out = fullfile(param.outdir, kernel_folders{i});
    CBIG_KRR_generate_kernels( feat_mat, param.sub_fold, kernel_out, param.ker_param)
end

%% step 3. Multi KRR
fprintf('# Step 3: Multi KRR.\n')
all_acc_metric = {'corr', 'COD', 'predictive_COD','MAE', 'MAE_norm',...
            'MSE','MSE_norm'};
for test_fold = 1: num_test_folds
    [acc_test,~,~, opt_stats] = CBIG_MultiKRR_FindOptTestAcc(param.outdir,...
    param.sub_fold,kernel_folders,param.outstem, param.num_inner_folds,test_fold,...
    param.domain,param.group_kernel, param.with_bias, param.threshold, param.ker_param, ...
    param.acc_metric);

    for ker_type = 1: length(param.ker_param)
        if test_fold == 1
            tmp =  acc_test(ker_type,:);
            acc{ker_type} = tmp;
            optimal_stats{ker_type} = opt_stats{ker_type};
        else
            tmp = acc_test(ker_type,:);
            acc{ker_type} = [acc{ker_type}; tmp];
            tmp_stats_struct = optimal_stats{ker_type};
            tmp_stats_struct_2 = opt_stats{ker_type};
            for kk = 1:length(all_acc_metric)
                tmp_stats = tmp_stats_struct(kk).value;
                tmp_stats_2 = tmp_stats_struct_2(kk).value;
                tmp_stats = [tmp_stats; tmp_stats_2];
                tmp_stats_struct(kk).value = tmp_stats;
            end
            optimal_stats{ker_type} = tmp_stats_struct;
               
        end
    end

end

save(fullfile(param.outdir, ['final_result_' param.outstem '.mat']), 'acc', 'optimal_stats')

    
