function [optimal_acc,optimal_stats] = CBIG_TRBPC_LRR_LpOCV_pick_optima(N_fold,N_score,data_dir, stem)

% [optimal_acc,optimal_stats] = CBIG_TRBPC_LRR_LpOCV_pick_optima(N_fold, data_dir, stem)
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
load(fullfile(data_dir,'behav_1','params','fold_1',['opt_results_' stem '.mat']));
stats_name = fieldnames(opt_stats);
N_stats = length(stats_name);

%% load all prediction results
optimal_acc = zeros(N_fold,N_score);
optimal_hyp = zeros(N_fold,N_score,2);
optimal_y_pred = cell(N_fold,N_score);
optimal_stats = opt_stats;
for i = 1:N_stats
    optimal_stats.(stats_name{i}) = zeros(N_fold,N_score);
end

for i = 1:N_fold
    for j = 1:N_score
        load(fullfile(data_dir,['behav_' num2str(j)],'params',['fold_' num2str(i)],['opt_results_' stem '.mat']));
        load(fullfile(data_dir,...
            ['behav_' num2str(j)],'params',['fold_' num2str(i)],['selected_parameters_' stem '.mat']));
        optimal_acc(i,j) = acc_corr;
        for k = 1:N_stats
            optimal_stats.(stats_name{k})(i,j) = opt_stats.(stats_name{k});
        end
        optimal_hyp(i,j,1) = curr_threshold;
        optimal_hyp(i,j,2) = curr_lambda;
        optimal_y_pred{i,j} = y_pred;
    end
end

save(fullfile(data_dir,['final_result_' stem '.mat']),'optimal_acc',...
    'optimal_hyp','optimal_y_pred','optimal_stats');