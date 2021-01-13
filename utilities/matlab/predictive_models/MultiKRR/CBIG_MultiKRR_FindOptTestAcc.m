function [opt_pred_acc,opt_y_pred,opt_hyp, opt_stats_cell] = CBIG_MultiKRR_FindOptTestAcc(data_dir,...
    sub_fold,kernel_folders,y_resid_stem, num_inner_folds,test_fold,...
    domain,group_kernel, with_bias, threshold, ker_param, acc_metric, gpso_dir)

% [acc_test,final_y_pred, opt_hyp, opt_stats_cell] = CBIG_MultiKRR_FindOptTestAcc(data_dir,...
%    sub_fold,kernel_folders,y_resid_stem, num_inner_folds,test_fold,...
%   domain,group_kernel, with_bias, threshold, ker_param, acc_metric, gpso_dir)
%
% This function performs multi KRR by using first using gaussian
% processes to search for the optimal hyperparameters before training the
% actual multi KRR model. This function assumes the following are stored in
% `data_dir`
%  1) Folder `y` which stores the regressed and original behavioral
%  measures. Generated using CBIG_crossvalid_regress_covariates_from_y.m 
%  2) Folders with name `kernel_folders{i}` whereby each of these folder
%  will store a kernel required to be trained by the multi KRR
%  model.Kernels in these folder are generated from
%  CBIG_KRR_generate_kernels
% 
% Inputs:
%   - data_dir 
%     A string. Contains the path to the directory storing the kernels and
%     regressed behaviors data
%
%   - sub_fold
%     The data split for cross-validation.
%     It is a num_test_folds x 1 structure with a field "fold_index".
%     sub_fold(i).fold_index is a #subjects x 1 binary vector, where 1
%     refers to the test subjects in the i-th fold; 0 refers to the
%     training subjects in the i-th fold.
%     If the user does not need cross-validation, the length of sub_fold
%     will be 1. In this case, sub_fold.fold_index has 3 unique values: 
%     0 - training set;
%     1 - test set;
%     2 - validation set.
%
%   - kernel_folders
%     Cell array of strings of size equal to the number of kernels to be passed in
%     to the multi-kernel regression model. Each element in the dictionary
%     contains the name of a sub folder in 'data_dir'. The different
%     sub_folders in 'data_dir' stores the kernel matrices to be used for
%     multi-kernel regression. For example, kernel_folders could be ofthe
%     form {'Kernel1', 'Kernel2'}. This means that the kernel matrices are
%     stored in [data_dir '/Kernel1'] and [data_dir '/Kernel2'].Each of the
%     kernel sub folder contains the kernel matrices generated from
%     CBIG_KRR_generate_kernels. For example, [data_dir
%     'Kernel1/FSM_innerloop'] and [data_dir 'Kernel1/FSM_test']. See the
%     descriptions in CBIG_KRR_generate_kernerls for more information
%     regarding the files stored in FSM_innerloop and FSM_test.
%     
%   - y_resid_stem
%     A string that was used to describe the .mat file of the target
%     variable y after regression. For example, if the filename is
%     <path_to_file>/y_regress_58behaviors.mat,
%     "y_resid_stem" will be '58behaviors'.
%     See the description of "outstem" parameter of function
%     CBIG_crossvalid_regress_covariates_from_y.m    
%
%   - num_inner_folds
%     A scalar, the number of inner-loop cross-validation folds for
%     hyperparameter selection.
%     If there are enough data so that cross-validation is not needed, set
%     num_inner_folds = 1.
%
%   - test_fold
%     A string or a scalar, current test fold.
%
%   - domain
%     A 1 x 2 vector. Stores the search domain of the hyperparameter
%     lambda.The first element of the vector stores the lower bound and the
%     second element stores the upper bound of the search domain.
%
%   - group_kernel (optional)
%     A 1 by #kernel_groups cell array. Each element in the cell array contains
%     a vector of size 1 x #kernels_in_the_group. For example if 4 kernels
%     are to be used to train the multi kernel model, but they can be
%     grouped such that 3 of the kernels uses the same hyperparameter
%     lambda for regularization, then I can pass in the 'group_kernel' as
%     {[1,3,4],[2]}. This means that kernels 1, 3 and 4 are regularized
%     using the same lambda while kernel 2 is regularized using a different
%     lambda. Do note that kernel 1 would correspond to the kernels stored
%     in the sub folder as specified by the first element of
%     'kernel_folders' so on and so forth.If group_kernel is empty ,i.e.
%     '', then it is assumed that each kernel will be regularised
%     independently.
%
%   - with_bias (optional)
%     A binary scalar(choose between 0 or 1).
%     - with_bias = 0 means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data. Default is
%     1.
%
%   - threshold (optional)
%     A scalar, it is a threshold required to determine the prediction (1
%     or 0) for binary target variables (e.g. sex). Ignore this parameter
%     if y is not binary. If y is binary but threhold is not passed in, the
%     default value is 0.5.
% 
%   - ker_param (optional)
%     A k x 1 structure with two fields: type and scale. k denotes the
%     number of kernel type.
%     ker_param(k).type is a string that specifies the type of kernel to be
%     used in the multi KRR model.
%                       Choose from
%                       'corr'        - Pearson's correlation;
%                       'Gaussian'    - Gaussian kernel;
%                       'Exponential' - exponential kernel.
%     ker_param(k).scale is a scalar specifying the scale of k-th kernel
%     (for Gaussian kernel or exponential kernel). If ker_param(k).type == 'corr',
%     ker_param(k).scale = NaN.
%     If "ker_param" is not passed in, only correlation kernel will be
%     calculated.Note that this function assumes that all the kernels used
%     to train a multi KRR model are of the same type. For example, if
%     ker_param(k).type is 'corr', and there are 2 kernel_folders, then
%     this function will train the multi KRR model using the correlation
%     kernels from each of the 2 kernel_folders.
%
%   - acc_metric (optional)
%     A string indicating the metric used to define accuracy that we wish to
%     optimise our hyperparameters on.
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
%   - gpso_dir (optional)
%     A string, the full path of the directory to store the cloned git
%     repository for Gaussian process code. Current script needs the
%     Gaussian-Process Surrogate Optimisation (GPSO;
%     https://github.com/jhadida/gpso) package to optimize the objective
%     function. This git repository and its dependency "deck"
%     (https://github.com/jhadida/deck) need to be both cloned into the
%     same folder, which is passed into current script through
%     params.gpso_dir.
%     Default is
%     $CBIG_CODE_DIR/external_packages/matlab/non_default_packages/Gaussian_Process
% 
% Outputs:
%   The outputs will be saved in [data_dir '/optimal_acc/fold_'
%   num2str(test_fold)]
%   - opt_pred_acc
%     A #k by #behavioral_measures matrix. This stores the optimal test
%     accuracies using correlation as metrics.
%
%   - opt_y_pred
%     A #k by #behavioral_measures cell array. Each element in the cell
%     array stores a #test_subject by 1 vector. This vector represents the
%     predicted behavioral measure for the given test subjects
%
%   - opt_hyp
%     A #k by #behavioral_measures cell array. Each element in the cell
%     array stores a 1 by #kernel_group vector. The first element in this
%     vector stores the optimal hyperlarameter (lambda) for the first
%     kernel_group and so on and so forth.
%
%   - opt_stats_cell
%     A #k sized cell array. Each element in the cell array stores a data
%     structure. There are 2 fields in the data structure, 'value' and
%     'description'. 'description' describes which accuracy metric the
%     'value' belongs to. 'value' will be a 1 x #behavioral_measures vector
%     storing the values of the optimal test prediction accuracies
%     corresponding to the specified acuracy metric. Within each field,
%     there will be a total number of 7 entries, each entry corresponding
%     to the possible acc_metrics.
%
% Example:
%     data_dir = '/data/users/xxx'
%     kernel_folders = {'kernel1','kernel2'};
%     y_resid_stem = 'Cognitive';
%     num_inner_folds = 20;
%     test_fold = 1;
%     domain = [0 20];
%     group_kernel = '';
%     
% 
%     [opt_pred_acc,opt_y_pred,opt_hyp] = CBIG_MultiKRR_FindOptTestAcc(data_dir,...
%     sub_fold,kernel_folders,y_resid_stem, num_inner_folds,test_fold,...
%     domain,group_kernel, with_bias, threshold, ker_param)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Yan Rui Tan and Ru(by) Kong

%% Set up
if(~exist('gpso_dir', 'var') || isempty(gpso_dir))
    gpso_dir = fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
        'non_default_packages', 'Gaussian_Process');
end
gpso_info = dir(fullfile(gpso_dir, 'deck'));
if(isempty(setdiff({gpso_info.name}, {'.', '..'})))
    curr_dir = pwd;
    cd(getenv('CBIG_CODE_DIR'))
    command = sprintf('git submodule update --init --recursive');
    system(command);
    cd(curr_dir)
end

addpath(fullfile(gpso_dir, 'gpso'))
addpath(fullfile(gpso_dir, 'deck'))
dk_startup();

%% Check vector Dimension
if size(domain,1) > 1
    error('The variable domain should contain a row vector')
end

%% Prep variables
num_kernel = length(kernel_folders);

if (isempty(group_kernel) || ~exist('group_kernel', 'var'))
    for i = 1:num_kernel
        group_kernel{i} = i;
    end
end

num_group_kernel =length(group_kernel);

if(ischar(test_fold))
    test_fold = str2num(test_fold);
end

if(ischar(num_inner_folds))
    num_inner_folds = str2num(num_inner_folds);
end

%% Load regressed behaviors. Contains 2 variables: y_resid and y_orig
if(~isempty(y_resid_stem))
    y_resid_stem = ['_' y_resid_stem];
end

y_resid_file = fullfile(data_dir, 'y', ['fold_' num2str(test_fold)], ...
['y_regress' y_resid_stem '.mat']);
load(y_resid_file)

%% set default hyperparamters if not passed in
if(~exist('ker_param', 'var') || strcmpi(ker_param, 'none'))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
end


if (isempty(acc_metric) || ~exist('acc_metric','var'))
    acc_metric = 'predictive_COD';
end


if (isempty(with_bias) || ~exist('with_bias','var'))
    with_bias = 1;
end
%% Check if behaviors are binary or continuous values
for i = 1:size(y_orig, 2)
    if(length(unique(y_orig(:,i))) > 2)
        % not binary
        bin_flag(i) = 0;
    else
        % binary case
        bin_flag(i) = 1;
    end
end
if(any(bin_flag==1) && any(bin_flag==0))
    error('Mixture of binary (e.g. sex) and continuous cases. Please run separately')
elseif(~any(bin_flag==1))    % all continuous
    bin_flag = 0;
    threshold = NaN;
else                         % all binary
    bin_flag = 1;
    if(~exist('threshold', 'var') || strcmpi(threshold, 'none'))
        threshold = 0.5;
    end
end

%% check if using cross-validation

 if(length(sub_fold)>1)
    % cross-validation
    num_valid_sub = [];
    kernel_innerloop_dir_suffix = fullfile('FSM_innerloop', ['fold_' num2str(test_fold)]);
    kernel_test_dir_suffix = fullfile('FSM_test', ['fold_' num2str(test_fold)]);
    sub_index_innerloop = sub_fold(test_fold).fold_index==0;
    sub_index_test = sub_fold(test_fold).fold_index==1;
else
    %Train validate test set
    assert(num_inner_loop==1,'"num_inner_loop" needs to be "1" as CV not performed')    
    num_valid_sub = length(find(sub_fold.fold_index==2));
    kernel_innerloop_dir_suffix = 'FSM_innerloop';
    kernel_test_dir_suffix = 'FSM_test';
    sub_index_innerloop = [find(sub_fold.fold_index==0); find(sub_fold.fold_index==2)];
    sub_index_test = sub_fold(test_fold).fold_index==1;
 end

%% Search for optimal hyperparam and find optimal test accuracies
outdir = fullfile(data_dir, 'optimal_acc', ['fold_' num2str(test_fold)]);
mkdir(fullfile(data_dir, 'optimal_acc'), ['fold_' num2str(test_fold)])

if(~exist(fullfile(outdir, ['acc' y_resid_stem '.mat']), 'file'))
    y_innloop_resid = y_resid(sub_index_innerloop, :);
    y_innloop_orig = y_orig(sub_index_innerloop, :);
    y_test_resid = y_resid(sub_index_test, :);
    y_test_orig = y_orig(sub_index_test, :);

    for k = 1:length(ker_param)
        fprintf('Kernel type: %s, scale: %f\n', ker_param(k).type, ker_param(k).scale)
        
        for i = 1:num_kernel
            if(strcmp(ker_param(k).type, 'corr'))
                ker = kernel_folders{i};
                kernel_innerloop = fullfile(data_dir, ker,...
                    kernel_innerloop_dir_suffix,...
                    ['FSM_' ker_param(k).type '.mat']);
                kernel_test = fullfile(data_dir, ker, kernel_test_dir_suffix,... 
                    ['FSM_' ker_param(k).type '.mat']);
            else
                kernel_innerloop = fullfile(data_dir, ker, kernel_innerloop_dir_suffix,...
                    ['FSM_' ker_param(k).type num2str(ker_param(k).scale) '.mat']);
                kernel_test = fullfile(data_dir, ker, kernel_test_dir_suffix,...
                    ['FSM_' ker_param(k).type num2str(ker_param(k).scale) '.mat']);
            end
            load(kernel_innerloop)
            FSM_innerloop(:,:,i) = FSM;
            clear FSM
            load(kernel_test)
            FSM_test_train(:,:,i) = FSM(sub_index_innerloop, sub_index_innerloop);
            FSM_test_pred(:,:,i) = FSM(sub_index_test, sub_index_innerloop);
            clear FSM
        end
    
        num_train_subject = size(FSM_innerloop,1);
        num_behav = size(y_innloop_resid,2);
        
        rng(1);
        cv_idx = cvpartition(num_train_subject,'kfold',num_inner_folds);
        
        domain = repmat(domain,num_group_kernel,1);
        
        for i = 1:num_behav
            % Use GSPO to search for optimal hyper parameter
            objfun = @(lambda) alpha_obj_function_multi(FSM_innerloop,lambda,...
                y_innloop_resid(:,i),y_innloop_orig(:,i), cv_idx,num_inner_folds,...
            with_bias, num_valid_sub, group_kernel, bin_flag, threshold, acc_metric);
            obj = GPSO();
            output =  obj.run(objfun, domain, 20,'tree',3);  
            
            % Extract out the optimal lambda and innerloop averaged
            % accuracies
            opt_lambda = output.sol.x;    
            acc_innloop(i) = output.sol.f;
            acc_train(i) = output.sol.f;
            
            % Prepare optimal hyperparameter to train the data  
            lambda_vect = zeros(1, num_kernel);
            for j = 1:length(group_kernel)
                kernels = group_kernel{j};
                lambda_vect(kernels) = opt_lambda(j);
            end
            
            % With the optimal hyperparameter, Use it to train and predict
            % behaviors on test set
            all_acc_metric = {'corr', 'COD', 'predictive_COD','MAE', 'MAE_norm',...
            'MSE','MSE_norm'};
            [~, opt_pred_acc(k,i),curr_opt_stats, opt_y_pred{k,i}] = CBIG_MultiKRR_TrainAndTest(...
                FSM_test_train, FSM_test_pred, y_innloop_resid(:,i), y_test_resid(:,i),...
                y_test_orig(:,i), lambda_vect, bin_flag, with_bias,threshold, acc_metric, 1);
            if i == 1
                opt_stats = curr_opt_stats;
            else
                for kk = 1:length(all_acc_metric)
                    tmp_opt_stats = opt_stats(kk).value;
                    tmp_opt_stats_2 = curr_opt_stats(kk).value;
                    tmp_opt_stats = [tmp_opt_stats, tmp_opt_stats_2];
                    opt_stats(kk).value = tmp_opt_stats;
                end
            end
            opt_stats_cell{k} = opt_stats;                
            opt_hyp{k,i} = lambda_vect;
            opt_acc_innloop{k,i} = acc_innloop;
            opt_acc_train{k,i} = acc_train;
            
        end
    end
    save(fullfile(outdir, ['acc' y_resid_stem '.mat']), 'opt_pred_acc','opt_hyp',...
        'opt_y_pred', 'opt_acc_innloop', 'opt_acc_train', 'opt_stats_cell')
else
    load(fullfile(outdir, ['acc' y_resid_stem '.mat']))
end

%% Remove Path
rmpath(fullfile(gpso_dir, 'gpso'))
rmpath(fullfile(gpso_dir, 'deck'))

end


function avg_acc = alpha_obj_function_multi(FSM_innloop,lambda,y_resid,...
    y_orig,cv_idx,num_innloop_fold,with_bias,num_valid_sub,group_kernel,...
    bin_flag,threshold, acc_metric)

% avg_acc = alpha_obj_function_multi(FSM_innloop,lambda,y_resid,...
%    y_orig,cv_idx,num_innloop_fold,with_bias,num_valid_sub,group_kernel,...
%    bin_flag,threshold, acc_metric)
%
% The Multi KRR uses Gaussian Processes to search for the optimal
% hyperparameter lambda. We are interested in finding the lambda that gives
% us the best accuracy (averaged across the innerloop folds).
% This function serves as the objective function to be passed into the GPML
% (Gaussian Process) module. This function takes in a specific 
% value of lambda and tabulates the averaged accuracy. By passing this
% function to the GPML module in matlab, the GPML uses this function to
% iteratively search for the best lambda that gives us the best accuracy.
% 
% Inputs:
%   - FSM_innloop 
%     A Nt x Nt x K matrix where Nt is number of training subjects and K is
%     number of kernels to be used (equivalent to the size of
%     kernel_folders in the main function). This matrix is the functional
%     similarity matrix used typically to be passed into KRR where each
%     element in the matrix stores the similarity of features between any 2
%     subjects.
%
%   - lambda
%     A 1 x #kernel_groups vector storing the lambda value for each kernel
%     group.
%
%   - y_resid
%     A Nt x 1 vector where Nt is number of training subjects. This vector
%     stores the target measure of all the triaining subjects after
%     regression of covariates.
%
%   - y_orig
%     A Nt x 1 vector where Nt is number of training subjects. This vector
%     stores the target measure of all the triaining subjects before
%     regression of covariates.
%
%   - num_innloop_fold
%     A scalar. Stores the number of innerloop folds you wish to have.
%
%   - with_bias
%     A binary scalar(choose between 0 or 1).
%     - with_bias = 0 means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
%
%   - num_valid_sub
%     A scalar or an empty variable. It should be a scalar when you are 
%     doing train,validation, test instead of cross validation. The scalar
%     represents the number of validation subjects you have. Otherwise,
%     input an empty matrix.
%
%   - group_kernel
%     A 1 by #kernel_groups cell array. Each element in the cell array contains
%     a vector of size 1 x #kernels_in_the_group. For example if 4 kernels
%     are to be used to train the multi kernel model, but they can be
%     grouped such that 3 of the kernels uses the same hyperparameter
%     lambda for regularization, then I can pass in the 'group_kernel' as
%     {[1,3,4],[2]}. This means that kernels 1, 3 and 4 are regularized
%     using the same lambda while kernel 2 is regularized using a different
%     lambda. Do note that kernel 1 would correspond to the kernels stored
%     in the sub folder as specified by the first element of
%     'kernel_folders' so on and so forth.If group_kernel is empty ,i.e.
%     '', then it is assumed that each kernel will be regularised
%     independently.
%
%   - bin_flag
%     A binary scalar. 1 means that the target measure are binary values. 0
%     means otherwise
%
%   - threshold
%     A scalar, it is a threshold required to determine the prediction (1
%     or 0) for binary target variables (e.g. sex). Ignore this parameter
%     if y is not binary. If y is binary but threhold is not passed in, the
%     default value is 0.5.
%
%   - acc_metric
%     A string indicating the metric used to define prediction accuracy. The
%     acc_metric prediction accuracy is the method used to calculate the 
%     prediction accuracy based on the calculated predicted target measures
%     value. 
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
% Outputs:
%   - avg_acc
%     A scalar. This is the averaged accuracy across all innerloop folds 
%     for a specific vector of lambdas.


for fold = 1:num_innloop_fold

    num_kernel = size(FSM_innloop,3);
    lambda_vect = zeros(1,num_kernel);
    for i = 1:length(group_kernel)
        kernels = group_kernel{i};
        lambda_vect(1,kernels) = lambda(i);
    end        
    
    if num_innloop_fold == 1
        idx_innloop_train = 1:(size(FSM_innloop,1) - num_valid_sub);
        idx_innloop_pred = (size(FSM_innloop,1) - num_valid_sub + 1):size(FSM_innloop,1);
    else
        idx_innloop_train = cv_idx.training(fold);
        idx_innloop_pred = cv_idx.test(fold);
    end

    FSM_innloop_train = FSM_innloop(idx_innloop_train,idx_innloop_train,:);
    FSM_innloop_pred = FSM_innloop(idx_innloop_pred,idx_innloop_train,:);

    y_resid_innloop_train = y_resid(idx_innloop_train,:);
    y_resid_innloop_pred = y_resid(idx_innloop_pred,:);
    y_true = y_orig(idx_innloop_pred,:);
  
    
    [innloop_acc(fold),~,~, ~] = CBIG_MultiKRR_TrainAndTest(FSM_innloop_train,...
        FSM_innloop_pred,y_resid_innloop_train,y_resid_innloop_pred,y_true,...
        lambda_vect,bin_flag,with_bias,threshold, acc_metric, 0);
    
    
end
avg_acc = mean(innloop_acc);

end
