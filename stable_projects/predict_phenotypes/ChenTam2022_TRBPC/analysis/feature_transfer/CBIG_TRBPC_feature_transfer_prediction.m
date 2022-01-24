function CBIG_TRBPC_feature_transfer_prediction(feature_path, pred_dir, PFM_average_dir, cluster, domain, outdir)

% CBIG_TRBPC_feature_transfer_prediction(feature_path, pred_dir, PFM_average_dir, cluster, domain, behav, outdir)
%
% This function prepares the files that will be used for the feature transfer
% learning in Chen & Tam 2021 paper. For each behavior, we compute the
% average PFM of all behaviors in a behavior domain excluding the behavior itself.
%
% Inputs:
%   - feature_path
%     Path of feature files that were used for prediction
%
%   - pred_dir
%     output directory of the multi-KRR prediction. We assume that files in
%     this directory are produced by CBIG_TRBPC_multiKRR_LpOCV_workflow.sh so
%     it follows certain folder structure.
%
%   - PFM_average_dir
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for hypothesis-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
%   - cluster
%     behavior clusters that was used to define the behavioral domain.
%     Choose from 'hypothesis' and 'datadriven'.
%
%   - domain
%     The behavioral domain. Choose from 'cog' for cognition, 'pers' for
%     personality and 'mh' for mental health.
%
%   - outdir
%     Output directory for the accuracy files
%
% Outputs:
%   Output files will be saved in outdir
%
% Written by Jianzhong Chen & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load(fullfile(PFM_average_dir,['PFM_domain_average_' cluster '_all_folds_' domain '.mat']), 'PFM_all_folds');
load(fullfile(pred_dir, 'no_relative_3_fold_sub_list.mat'), 'sub_fold');

% default percentage of choosing edges is 10
perc = 10;
feat_mat = 0;
load(feature_path,'feature_files');
for i = 1:length(feature_files)
    curr_feature = load(feature_files{i});
    names = fieldnames(curr_feature);
    feat_mat = feat_mat + curr_feature.(names{1});
end

N_edge = 419*418/2;
acc_score_fold = zeros(36,120);

for i = 1:120
    load(fullfile(pred_dir, 'y', ['fold_' num2str(i)], 'y_regress_all_score.mat'));
    for behav = 1:36
        test_ind = sub_fold(i).fold_index;
        importance_curr_task = 0;
        for k = 1:4
            importance_curr_task = importance_curr_task + PFM_all_folds((k-1)*N_edge+1:k*N_edge,i,behav);
            
        end
        [~,ind] = sort(abs(importance_curr_task),'Descend');
        ind_top = ind(1:floor(length(importance_curr_task)*perc/100));
        selected_FC = feat_mat(ind_top,:);
        pos_ind = importance_curr_task(ind_top) > 0;
        neg_ind = importance_curr_task(ind_top) < 0;
        feature = sum(selected_FC(pos_ind,:),1) - sum(selected_FC(neg_ind,:),1);
        y_test = y_resid(test_ind,behav);
        x_test = feature(test_ind)';
        acc_score_fold(behav,i) = corr(x_test,y_test);
    end
end
save(fullfile(outdir, ['feature_transfer_acc_' cluster '_' domain '_' num2str(perc) 'perc.mat']), 'acc_score_fold');
end
