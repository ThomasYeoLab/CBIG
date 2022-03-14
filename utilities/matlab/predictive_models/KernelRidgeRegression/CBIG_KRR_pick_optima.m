function [optimal_acc, optimal_kernel, optimal_lambda, optimal_threshold] = ...
    CBIG_KRR_pick_optima( sub_fold, data_dir, stem, bin_flag, ...
    ker_param, lambda_set, threshold_set )

% [optimal_acc, optimal_kernel, optimal_lambda, optimal_threshold] = ...
%     CBIG_KRR_pick_optima( sub_fold, data_dir, stem, bin_flag, ...
%     ker_param, lambda_set, threshold_set )
% 
% Based on the inner-loop cross-validation accuracies (or validation
% accuracies), this function selects the optimal hyperparameters: kernel
% parameters, regularization parameter lambda, and/or the threshold to
% calculate accuracies (for binary target variables). After the optimal
% hyperparameters are selected, the corresponding test accuracies will be
% read out and saved in a final results (.mat) file.
% 
% Inputs:
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
%     Full path of the input/output directory. Assumptions of subfolders:
%     (1) [data_dir '/innerloop_cv/fold_' num2str(test_fold)]: stores
%         inner-loop cross-validation results for each test fold. See
%         CBIG_KRR_innerloop_cv_allparams.m
%     (2) [data_dir '/test_cv/fold_' num2str(test_fold)]: stores
%         training-test cross-validation results for each test fold. See
%         CBIG_KRR_test_cv_allparams.m
%     A .mat file [data_dir '/final_result_' stem '.mat'] will be created
%     to save the optimal hyperparameters and corresponding accuracy and
%     predicted scores.
% 
%   - stem
%     A string appended to specify the input files. For example, if the
%     inner-loop CV and test CV accuracy filename, as output from
%     CBIG_KRR_innerloop_cv_allparam.m, has the format
%     <path_to-file>/acc_58behaviors.mat, then stem = '58behaviors'.
%     If accuracy input files are <path_to-file>/acc.mat, then stem = ''.
% 
%   - bin_flag
%     0 or 1. 1 means the variables to be predicted (i.e. y_orig) are
%     binary. 0 means none of them are binary.
%     If your target variables contain both binary and non-binary measures,
%     you need to deal with they separately with different bin_flag.
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
%     affect the accuracy of the target variable (binary) predictions. To be
%     used when the target variable to predict is binary (eg. sex). The
%     threshold is used as a point to divide the prediction into the two
%     classes (eg. Male or Female). "Thresholds" should be between -1 and
%     1. 
%     If "threshold_set" is not passed in or is 'NONE', it will be set as
%     default: [-1:0.1:1].
% 
% Outputs:
%   - optimal_acc
%     A #TestFolds x #TargetVariable matrix of optimal accuracy of each
%     test fold and each target measure predicted.
%  
%   - optimal_kernel
%     A structure with length of #TestFolds.
%     optimal_kernel(i).type denotes the optimal kernel's type
%     ('corr', 'Gaussian', or 'exponential') of i-th test fold.
%     optimal_kernel(i).scale denotes the optimal kernel's scale (a scalar
%     when optimal_kernel(i).type == 'Gaussian' or 'exponential'; NaN when
%     optimal_kernel(i).type == 'corr').
% 
%   - optimal_lambda
%     A #TestFolds x #TargetVariable vector of optimal lambda of each
%     test fold and each measure to predict.
% 
%   - optimal_threshold
%     A #TestFolds x #TargetVariable vector of optimal threshold for the
%     case of bin_flag == 1 of each test fold and each measure to predict.
%     If bin_flag == 0, then every entry of optimal_threshold is NaN.
%
%   - y_predict_concat (saved in final result mat file)
%     A #TotalSubjects x #TargetVariable matrix of the predicted target variable for 
%     each subject.
%
%   - optimal_stats (saved in final result mat file)
%     A struct containing a #TestFolds x #TargetVariable matrix containing the 
%     prediction for each possible accuracy metric (eg. corr, MAE, etc).
%
%   - y_pred_train (saved in final result mat file)
%     A cell array of size equal to #TestFolds. Each cell contains a matrix of size
%     #TrainingSubjects x #TargetVariable element of the predicted target variable 
%     for the training subjects of each fold.
%     
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Jingwei Li and Ru(by) Kong

%% set default hyperparamters if not passed in
if(~exist('ker_param', 'var') || strcmpi(ker_param, 'none'))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
end

if(~exist('lambda_set', 'var') || strcmpi(lambda_set, 'none'))
    lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 ...
        5 10 15 20];
end

if(bin_flag==1)
    if(~exist('threshold_set', 'var') || strcmpi(threshold_set, 'none') || isempty(threshold_set))
        threshold_set = [-1:0.1:1];
    end
else
    threshold_set = NaN;
end

if(~isempty(stem))
    stem = ['_' stem];
end
num_test_folds = length(sub_fold);
metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};

%%
y_predict_concat = [];
for test_fold = 1:num_test_folds
    innerloop = fullfile(data_dir, 'innerloop_cv', ['fold_' num2str(test_fold) '/acc' stem '.mat']);
    testloop = fullfile(data_dir, 'test_cv', ['fold_' num2str(test_fold) '/acc' stem '.mat']);
    innerloop = load(innerloop);
    testloop = load(testloop);
    
    innerloop_loss = CBIG_KRR_reorganize_acc(innerloop.loss);
    test_acc = CBIG_KRR_reorganize_acc(testloop.acc);
    
    for i = 1:size(innerloop_loss, 4)
        curr_innerloop_loss = innerloop_loss(:,:,:,i);
        curr_test_acc = test_acc(:,:,:,i);
        
        [~,I] = min(curr_innerloop_loss(:));
        [FSM_idx, lambda_idx, threshold_idx]=ind2sub(size(curr_innerloop_loss),I);
        
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
        y_pred_train{test_fold}(:,i) = testloop.y_pred_train{FSM_idx, lambda_idx, threshold_idx}{i};
        if(num_test_folds==1)
            y_predict_concat(:,i) = testloop.y_p{FSM_idx, lambda_idx, threshold_idx}{i};
        else
            y_predict_concat(sub_fold(test_fold).fold_index==1,i) = ...
                testloop.y_p{FSM_idx, lambda_idx, threshold_idx}{i};
        end
    end
end

save(fullfile(data_dir, ['final_result' stem '.mat']), 'y_predict_concat', 'optimal_kernel', ...
    'optimal_lambda', 'optimal_threshold', 'optimal_acc', 'optimal_stats', 'y_pred_train')


end


function acc_out = CBIG_KRR_reorganize_acc(acc)

% CBIG_KRR_reorganize_acc(acc)
% acc: 
%   A cell array with dimension of #kernels x #lambda x #threshold. Each
%   array is a 1 x #VariableToPredict matrix
% acc_out:
%   A #kernels x #lambda x #threshold x #VariableToPredict matrix.

for k = 1:size(acc, 1)
    for l = 1:size(acc, 2)
        for t = 1:size(acc, 3)
            curr_acc = acc{k,l,t};
            curr_acc = reshape(curr_acc, [1 1 1 length(curr_acc)]);
            acc_out(k,l,t,:) = curr_acc;
        end
    end
end

end

