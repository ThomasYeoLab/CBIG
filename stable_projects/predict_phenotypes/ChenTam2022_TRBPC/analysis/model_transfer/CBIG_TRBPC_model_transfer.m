function CBIG_TRBPC_model_transfer(pred_dir, hypothesis_ind, datadriven_ind, outdir)

% CBIG_TRBPC_model_transfer(pred_dir, cluster, outdir)
%
% This function prepares the files that will be used for the feature transfer
% learning in Chen & Tam 2021 paper. For each behavior, we compute the
% average PFM of all behaviors in a behavior domain excluding the behavior itself.
%
% Inputs:
%   - pred_dir
%     output directory of the multi-KRR prediction. We assume that files in
%     this directory are produced by CBIG_TRBPC_multiKRR_LpOCV_workflow.sh so
%     it follows certain folder structure.
%
%   - hypothesis_ind
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for hypothesis-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
%   - datadriven_ind
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for data-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
%   - outdir
%     Output directory for the accuracy files
%
% Outputs:
%   Output files will be saved in outdir
%
% Written by Jianzhong Chen & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load(fullfile(pred_dir, 'final_result_all_score.mat'));
load(fullfile(pred_dir, 'no_relative_3_fold_sub_list.mat'), 'sub_fold');
models = {'cog','pers','mh'};

%% for hypothesis-driven behavioral clusters
acc_score_model_fold = zeros(36,3,120);
for d = 1:length(models)
    ind_curr_domain = hypothesis_ind.(models{d});
    acc_all = zeros(length(ind_curr_domain),length(sub_fold));
    for i = 1:length(sub_fold)
        load(fullfile(pred_dir, 'y', ['fold_' num2str(i)], 'y_regress_all_score.mat'));
        y_pred = optimal_y_pred{i};
        for j = 1:36
            y = y_resid(sub_fold(i).fold_index,j);
            transfer_pred = zeros(size(y));
            for k = 1:length(ind_curr_domain)
                if ind_curr_domain(k) ~= j
                    transfer_pred = y_pred{ind_curr_domain(k)} + transfer_pred;
                end
            end
            transfer_pred = transfer_pred/sum(ind_curr_domain ~= j);
            acc_all(j,i) = corr(transfer_pred,y);
        end
    end
    acc_score_model_fold(:,d,:) = acc_all;
end
save(fullfile(outdir, 'acc_score_model_fold_hypothesis'), 'acc_score_model_fold');

%% for data-driven behavioral clusters
acc_score_model_fold = zeros(36,3,120);
for d = 1:length(models)
    ind_curr_domain = datadriven_ind.(models{d});
    acc_all = zeros(length(ind_curr_domain),length(sub_fold));
    for i = 1:length(sub_fold)
        load(fullfile(pred_dir, 'y', ['fold_' num2str(i)], 'y_regress_all_score.mat'));
        y_pred = optimal_y_pred{i};
        for j = 1:36
            y = y_resid(sub_fold(i).fold_index,j);
            transfer_pred = zeros(size(y));
            for k = 1:length(ind_curr_domain)
                if ind_curr_domain(k) ~= j
                    transfer_pred = y_pred{ind_curr_domain(k)} + transfer_pred;
                end
            end
            transfer_pred = transfer_pred/sum(ind_curr_domain ~= j);
            acc_all(j,i) = corr(transfer_pred,y);
        end
    end
    acc_score_model_fold(:,d,:) = acc_all;
end
save(fullfile(outdir, 'acc_score_model_fold_datadriven'), 'acc_score_model_fold');
