function [acc, pred_stats] = CBIG_KRR_test_cv_allparams_LITE( test_fold, sub_fold, ...
    data_dir, y_resid_stem, with_bias, kernel, ker_param, lambda_set, threshold_set )

% [acc] = CBIG_KRR_test_cv_allparams_LITE( test_fold, sub_fold, ...
%     data_dir, y_resid_stem, ker_param, lambda_set, threshold_set )
% 
% For a given test fold, this function calculates the test accuracies of
% all target variables with all possible combinations of hyperparameters,
% including the kernel parameters, the regularizaiton parameter lambda,
% and/or the threshold used for binary target variable.
% 
% This is a lite version of CBIG_KRR_test_cv_allparams.m. In this version,
% the kernels are not computed and saved separately for each inner-loop
% fold and each outer training-test fold. Instead, the kernel is only
% computed once across all the subjects, and for each fold, the scripts
% will grab the corresponding entries on the fly.
% 
% Note: if your target variables contain both binary and non-binary
% measures, please deal with them separately.
% 
% Inputs:
%   - test_fold
%     A string or a scalar, current test fold.
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
%   - data_dir
%     Full path of the input/output directory.
%     A subfolder [data_dir '/test_cv/fold_' num2str(test_fold)] will
%     be created to save the test cross-validation results.
% 
%   - y_resid_stem
%     A string that was used to describe the .mat file of the target
%     variable y after regression. For example, if the filename is
%     <path_to_file>/y_regress_58behaviors.mat,
%     "y_resid_stem" will be '58behaviors'.
%     See the description of "outstem" parameter of function
%     CBIG_crossvalid_regress_covariates_from_y.m
%   
%   - with_bias
%     A (choose from 0 or 1).
%     - with_bias = 0 means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
% 
%   - kernel
%     a #AllSubjects x #AllSubjects kernel matrix (The ordering of subjects 
%     follows the original subject list, i.e. the ordering in "feature_mat"
%     or "similarity_mat".)
%     A 'none' will be passed in, if it was saved in [outdir '/FSM']. 
%     Algorithm will load it from the file.
%
%   - ker_param (optional)
%     A K x 1 structure with two fields: type and scale. K denotes the
%     number of kernels.
%     ker_param(k).type is a string that specifies the type of k-th kernel.
%                       Choose from
%                       'corr'        - Pearson's correlation;
%                       'Gaussian'    - Gaussian kernel;
%                       'Exponential' - exponential kernel.
%     ker_param(k).scale is a scalar specifying the scale of k-th kernel
%     (for Gaussian kernel or exponential kernel). If ker_param(k).type == 'corr',
%     ker_param(k).scale = NaN.
%     If "ker_param" is not passed in, only correlation kernel will be
%     calculated.
% 
%   - lambda_set (optional)
%     A vector of numbers for grid search of lambda (the regularization
%     parameter). If "lambda_set" is not passed in or is 'NONE', it will be
%     set as default:
%     [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
%        5 10 15 20]
% 
%   - threshold_set (optional)
%     A vector of numbers (thresholds). Different thresholds used will
%     affect the accuracy of the target variable (binary) predictions. To
%     be used when the target variable to predict is binary (e.g. sex). The
%     threshold is used as a point (between -1 and 1) to divide the
%     prediction into the two classes (e.g. Male or Female).
%     If "threshold_set" is not passed in or is 'NONE', it will be set as
%     default: [-1:0.1:1].
% 
% Outputs:
%   - acc
%     A cell with dimension of #kernels x #lambda x #thresholds of
%     test CV accuracy for each hyperparameter combination.
%
%   - pred_stats
%     A cell array with dimension of #kernels x #lambda x #thresholds. Each
%     entry of the cell array is an matrix storing the prediction statistics
%     defined by 'corr', 'COD', 'predictive_COD', 'MSE', 'MSE_norm',
%     'MAE', 'MAE_norm'
%
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% setting up
if(ischar(test_fold))
    test_fold = str2num(test_fold);
end

if(~isempty(y_resid_stem))
    y_resid_stem = ['_' y_resid_stem];
end
y_resid_file = fullfile(data_dir, 'y', ['fold_' num2str(test_fold)], ...
    ['y_regress' y_resid_stem '.mat']);
load(y_resid_file)

% set default hyperparamters if not passed in
if(~exist('ker_param', 'var') || strcmpi(ker_param, 'none'))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
end

if(~exist('lambda_set', 'var') || strcmpi(lambda_set, 'none'))
    lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
        5 10 15 20];
end

% set the interested statistics
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

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
    error('Mixture of binary (e.g. sex) and continuous cases. Please run them separately.')
elseif(~any(bin_flag==1))    % all y are continuous
    bin_flag = 0;
    threshold_set = NaN;
else                         % all y are binary
    bin_flag = 1;
    if(~exist('threshold_set', 'var') || strcmpi(threshold_set, 'none') || isempty(threshold_set))
        threshold_set = [-1:0.1:1];
    end
end

kernel_dir = fullfile(data_dir, 'FSM');
% check if using cross-validation
if(length(sub_fold)>1)
    % cross-validation
    train_ind = sub_fold(test_fold).fold_index==0;
    test_ind = sub_fold(test_fold).fold_index==1;
    train_ind_y = train_ind;
    test_ind_y = test_ind;
else
    train_ind = find(sub_fold.fold_index==0);
    test_ind = find(sub_fold.fold_index==1);
    train_ind_y = train_ind;
    test_ind_y = test_ind;
end


%% test CV for each parameter combination
outdir = fullfile(data_dir, 'test_cv', ['fold_' num2str(test_fold)]);

if(strcmpi(kernel, 'none'))
    clear kernel
    for k = 1:length(ker_param)
        if(strcmp(ker_param(k).type, 'corr'))
            kfile = [kernel_dir '/FSM_' ker_param(k).type '.mat'];
        else
            kfile = [kernel_dir '/FSM_' ker_param(k).type '_' num2str(ker_param(k).scale) '.mat'];
        end
        kernel{k} = load(kfile);
    end
end

if(~exist(fullfile(outdir, ['acc' y_resid_stem '.mat']), 'file'))
    curr_y_resid_train = y_resid(train_ind_y, :);
    curr_y_resid_test = y_resid(test_ind_y, :);
    curr_y_orig_test = y_orig(test_ind_y, :);
    
    for k = 1:length(ker_param)
        fprintf('Kernel type: %s, scale: %f\n', ker_param(k).type, ker_param(k).scale)
        
        FSM_train = kernel{k}.FSM(train_ind, train_ind);
        FSM_test = kernel{k}.FSM(test_ind, train_ind);
        
        for l = 1:length(lambda_set)
            lambda = lambda_set(l);
            fprintf('  lambda: %f\n', lambda);
            
            for t = 1:length(threshold_set)
                threshold = threshold_set(t);
                fprintf('    threshold: %f\n', threshold)
                [y_p{k,l,t}, y_t{k,l,t}, acc{k,l,t}, pred_stats{k,l,t},y_pred_train{k,l,t}] = ...
                    CBIG_KRR_test_cv( bin_flag,FSM_train, FSM_test, curr_y_resid_train, ...
                    curr_y_resid_test, curr_y_orig_test, with_bias, lambda, threshold, metrics );
                clear threshold
            end
            clear lambda
        end
        clear FSM_train FSM_test
    end
    
    if(~exist(outdir, 'dir'))
        mkdir(outdir)
    end
    save(fullfile(outdir, ['acc' y_resid_stem '.mat']), 'acc', 'y_p', 'y_t', 'pred_stats', 'y_pred_train');
else
    load(fullfile(outdir, ['acc' y_resid_stem '.mat']))
    if(size(acc,1)~=length(ker_param) || size(acc,2)~=length(lambda_set) || size(acc,3)~=length(threshold_set))
        error('Wrong size of ''acc'' in existing file %s. Please remove it and rerun this script.', ...
            fullfile(outdir, ['acc' y_resid_stem '.mat']))
    end
end


end

