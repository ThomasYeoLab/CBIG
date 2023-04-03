function CBIG_GradPar_KRR_optimize_res_wrapper(params)

% CBIG_GradPar_KRR_optimize_res_wrapper(params)
%
% This script optimizes the resolution for KRR prediction using multiple resolutions for generating FC matrices.
% This script can also be used to optimize other parameter for KRR prediction. Specifically, given a set of KRR
% predictions of different approaches using the same set of lambda values, this script can be used to find the
% optimal result for each behavior across all approaches. In Kong2023, we generated a set of prediction results
% using different resolutions of FC matrices. This script is used to find the optimal resolution for each behavior.
%
% Please note that the prediction parameters such as num_splits, num_folds, num_behaviors, kernel parameters,
% and hyperparameter range should be consistent across the projects under project_set.
%
% Input:
%   - params: a structure with the following fields
%     - params.project_set: 
%       a cell array of strings. Each string is the path to the project which contains the esults of KRR prediction.
%       For example, if you want to optimize the resolution for Kong2021, then params.project_set = 
%       {'<path>/Kong2021/100', '<path>/Kong2021/200','<path>/Kong2021/300'}, where there are 3 resolutions for
%       Kong2021: 100, 200, and 300. The KRR prediction results for each resolution are stored in subfolders
%       under folder Kong2021. If the KRR prediction results were stored in the following folders: Kong2021_100,
%       Kong2021_200, and Kong2021_300, then the input can be params.project_set = {'<path>/Kong2021_100',
%       '<path>/Kong2021_200','Kong2021_300'}.
%
%     - params.nsplit:
%       a string. The current index of random splits of the cross-validation. For example, if you perform 100 repeats of
%       cross-validation, the current one is the first split, then params.nsplit = '1'. If there is no split, then
%       params.nsplit = ''.
%
%     - params.num_folds:
%       a string. The number of folds of the outer-loop cross-validation. For example, if you perform 20-fold
%       cross-validation, then params.num_folds = '20'.
%
%     - params.num_behaviors:
%       a string. The number of behaviors. For example, if you have 61 behaviors to be predicted, then
%       params.num_behaviors = '61'.
%
%     - params.bin_flag: (optional)
%       0 or 1. 1 means the variables to be predicted (i.e. y_orig) are binary. 0 means none of them are binary.
%       If your target variables contain both binary and non-binary measures, you need to deal with they
%       separately with different bin_flag. By default, bin_flag is 0.
%
%     - params.ker_param: (optional)
%       A Kx1 structure with two fields: type and scale. K is the number of kernels.
%       ker_param(k).type is a string that specifies the type of k-th kernel.
%       Choose from
%       'corr'        - Pearson's correlation;
%       'Gaussian'    - Gaussian kernel;
%       'Exponential' - exponential kernel.
%       ker_param(k).scale is a scalar specifying the scale of kernel
%       (for Gaussian kernel or exponential kernel). If ker_param(k).type == 'corr',
%       ker_param(k).scale = NaN.
%       If "ker_param" is not passed in, only correlation kernel will be
%       calculated.
%
%     - params.lambda_set: (optional)
%       A vector of numbers for grid search of lambda (the regularization parameter). If "lambda_set"
%       is not passed in or is 'NONE', it will be set as default:
%       [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
%        5 10 15 20]
% 
%     - params.threshold_set: (optional)
%       A vector of numbers (thresholds). Different thresholds used will affect the accuracy of the target
%       variable (binary) predictions. To be used when the target variable to predict is binary (eg. sex). The
%       threshold is used as a point to divide the prediction into the two classes (eg. Male or Female).
%       "Thresholds" should be between -1 and 1. If "threshold_set" is not passed in or is 'NONE', it will be
%       set as default: [-1:0.1:1].
%
%     - params.outstem:
%       a string. The name of the output file name. To keep the structure to be the same as the KRR prediction
%       results of the project_set, the output file name should be consistent with project_set. For example, if
%       the outstem of Kong2021/100 is '58_behaviors_3_components', then you should also use this string as the
%       params.outstem.
%
%     - params.out_dir:
%       a string. The path to the output folder. For example, if you want to optimize the resolution for Kong2021,
%       then you can give params.outname = '<path>/Kong2021_opt_res'.
%
% Output:
%   The output results will be stored in the folder params.out_dir. The output results will be saved in the same
%   structure as the KRR prediction results of the project_set.
%     - optimal_acc
%       A #TestFolds x #TargetVariable matrix of optimal accuracy of each test fold and each target measure
%       predicted.
%
%     - optimal_kernel
%       A structure with length of #TestFolds.
%       optimal_kernel(i).type denotes the optimal kernel's type ('corr', 'Gaussian', or 'exponential') of i-th
%       test fold.
%       optimal_kernel(i).scale denotes the optimal kernel's scale (a scalar when optimal_kernel(i).type == 
%       'Gaussian' or 'exponential'; NaN when
%       optimal_kernel(i).type == 'corr').
% 
%     - optimal_lambda
%       A #TestFolds x #TargetVariable vector of optimal lambda of each test fold and each measure to predict.
% 
%     - optimal_threshold
%       A #TestFolds x #TargetVariable vector of optimal threshold for the case of bin_flag == 1 of each test
%       fold and each measure to predict. If bin_flag == 0, then every entry of optimal_threshold is NaN.
%
%     - optimal_project
%       A #TestFolds x #TargetVariable cell of optimal project name of each test fold and each target variable.
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% set default hyperparamters if not passed in
% Check the struct params contains a field ker_param:
if(~isfield(params, 'ker_param') || isempty(params.ker_param))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
else
    ker_param = params.ker_param;
end
% Check the struct params contains a field lambda_set:
if(~isfield(params, 'lambda_set') || isempty(params.lambda_set) || strcmpi(params.lambda_set, 'none'))
    lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
        5 10 15 20];
else
    lambda_set = params.lambda_set;
end
% Check the struct params contains a field bin_flag:
if(~isfield(params, 'bin_flag') || isempty(params.bin_flag))
    bin_flag = 0;
else
    bin_flag = params.bin_flag;
end
% Check the struct params contains a field threshold_set:
if(bin_flag==1)
    if(~isfield(params, 'threshold_set') || isempty(params.threshold_set) || strcmpi(params.threshold_set, 'none'))
        threshold_set = [-1:0.1:1];
    else
        threshold_set = params.threshold_set;
    end
else
    threshold_set = NaN;
end
% Check the struct params contains a field outstem:
if(~isfield(params, 'outstem') || isempty(params.outstem))
    stem = '_';
else
    stem = ['_' params.outstem];
end
num_test_folds = str2num(params.num_folds);
sub_fold = params.sub_fold;

metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
project_set = params.project_set;
nsplit = params.nsplit;

y_predict_concat = [];
for test_fold = 1:num_test_folds
    for ns = 1:length(project_set)
        out_dir = fullfile(project_set{ns}, num2str(nsplit)); 
        innerloop = fullfile(out_dir, 'innerloop_cv', ['fold_' num2str(test_fold)], ['acc' stem '.mat']);
        testloop = fullfile(out_dir, 'test_cv', ['fold_' num2str(test_fold)], ['acc' stem '.mat']);
        innerloop = load(innerloop);
        testloop = load(testloop);
        innerloop_loss = CBIG_KRR_reorganize_data(innerloop.loss);
        test_acc = CBIG_KRR_reorganize_data(testloop.acc);
        for i = 1:size(innerloop_loss, 4)
            curr_innerloop_loss = innerloop_loss(:,:,:,i);
            curr_test_acc = test_acc(:,:,:,i);
            [min_val,I] = min(curr_innerloop_loss(:));
            min_val_projects(ns,i) = min_val;

        end
    end
    [~, min_project] = min(min_val_projects);

    for i = 1:size(innerloop_loss, 4)
        out_dir = fullfile(project_set{min_project(i)},num2str(nsplit)); 
        innerloop = fullfile(out_dir, 'innerloop_cv', ['fold_' num2str(test_fold)],['acc' stem '.mat']);
        testloop = fullfile(out_dir, 'test_cv', ['fold_' num2str(test_fold)], ['acc' stem '.mat']);
        innerloop = load(innerloop);
        testloop = load(testloop);
        innerloop_loss = CBIG_KRR_reorganize_data(innerloop.loss);
        test_acc = CBIG_KRR_reorganize_data(testloop.acc);

        curr_innerloop_loss = innerloop_loss(:,:,:,i);
        curr_test_acc = test_acc(:,:,:,i);
        
        [~,I] = min(curr_innerloop_loss(:));
        [FSM_idx, lambda_idx, threshold_idx]=ind2sub(size(curr_innerloop_loss),I);

        optimal_project{test_fold, i} = project_set{min_project(i)};
        optimal_kernel(test_fold, i).type = ker_param(FSM_idx).type;
        optimal_kernel(test_fold, i).scale = ker_param(FSM_idx).scale;
        optimal_lambda(test_fold, i) = lambda_set(lambda_idx);
        optimal_threshold(test_fold, i) = threshold_set(threshold_idx);
        optimal_acc(test_fold, i) = curr_test_acc(FSM_idx, lambda_idx, threshold_idx);
        for metric_ind = 1:length(metrics)
            optimal_stats.(metrics{metric_ind})(test_fold,i) = testloop.pred_stats{FSM_idx,...
                lambda_idx, threshold_idx}(metric_ind,i);
        end
        y_predict{test_fold}(:,i) = testloop.y_p{FSM_idx, lambda_idx, threshold_idx}{i};
        if(num_test_folds==1)
            y_predict_concat(:,i) = testloop.y_p{FSM_idx, lambda_idx, threshold_idx}{i};
        else
            y_predict_concat(sub_fold(test_fold).fold_index==1,i) = ...
                testloop.y_p{FSM_idx, lambda_idx, threshold_idx}{i};
        end
    end
    
end
fin_dir = fullfile(params.out_dir, num2str(nsplit));
mkdir(fin_dir);
save(fullfile(fin_dir, ['final_result' stem '.mat']), 'y_predict_concat', 'optimal_kernel', ...
    'optimal_lambda', 'optimal_threshold', 'optimal_acc', 'optimal_stats','optimal_project')


end


function data_out = CBIG_KRR_reorganize_data(data)

% CBIG_KRR_reorganize_data(data)
% data: 
%   A cell array with dimension of #kernels x #lambda x #threshold. Each
%   array is a 1 x #VariableToPredict matrix
% data_out:
%   A #kernels x #lambda x #threshold x #VariableToPredict matrix.

for k = 1:size(data, 1)
    for l = 1:size(data, 2)
        for t = 1:size(data, 3)
            curr_data = data{k,l,t};
            curr_data = reshape(curr_data, [1 1 1 length(curr_data)]);
            data_out(k,l,t,:) = curr_data;
        end
    end
end

end