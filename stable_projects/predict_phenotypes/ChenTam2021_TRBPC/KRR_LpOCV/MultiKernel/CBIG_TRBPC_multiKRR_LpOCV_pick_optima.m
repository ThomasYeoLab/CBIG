function [optimal_acc,optimal_stats] = CBIG_TRBPC_multiKRR_LpOCV_pick_optima(N_fold, data_dir, stem)

% [optimal_acc,optimal_stats] = CBIG_TRBPC_multiKRR_LpOCV_pick_optima(N_fold, data_dir, stem)
%
% Gather the prediction results of all cross-validation folds.
%
% Inputs:
%   - N_fold
%     Scalar. Number of cross-validation folds
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
% Outputs:
%   - optimal_acc
%     A #TestFolds x #TargetVariable matrix of optimal accuracy (as correlation) of each
%     test fold and each target measure predicted.
%
%   - optimal_stats
%     A struct. Each field is a A #TestFolds x #TargetVariable matrix of
%     prediction stats
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% get data dimension
load(fullfile(data_dir,'optimal_acc','fold_1',['acc_' stem '.mat']));
opt_stats = opt_stats_cell{1};
N_score = length(opt_pred_acc);
N_stats = length(opt_stats_cell{1});

%% load all prediction results
optimal_acc = zeros(N_fold,N_score);
optimal_hyp = cell(N_fold,1);
optimal_y_pred = cell(N_fold,1);
stats_name = cell(N_stats,1);
for i = 1:N_stats
    stats_name{i} = opt_stats(i).description;
    optimal_stats.(stats_name{i}) = zeros(N_fold,N_score);
end
    
for i = 1:N_fold
    load(fullfile(data_dir,'optimal_acc',['fold_' num2str(i)],['acc_' stem '.mat']));
    optimal_acc(i,:) = opt_pred_acc;
    opt_stats = opt_stats_cell{1};
    for j = 1:N_stats
        optimal_stats.(stats_name{j})(i,:) = opt_stats(j).value;
    end
    optimal_hyp{i} = opt_hyp;
    optimal_y_pred{i} = opt_y_pred;
end

save(fullfile(data_dir,['final_result_' stem '.mat']),'optimal_acc',...
    'optimal_hyp','optimal_y_pred','optimal_stats');