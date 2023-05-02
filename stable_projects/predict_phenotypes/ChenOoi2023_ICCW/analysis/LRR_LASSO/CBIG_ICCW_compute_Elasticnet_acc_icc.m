function CBIG_ICCW_compute_Elasticnet_acc_icc(model, result_dir, output_dir)

% function CBIG_ICCW_compute_Elasticnet_acc_icc(model, result_dir, output_dir)
%
% This function collates the accuracy results, and computes the ICC for the Haufe
% transformed weights, as well as the original regression weights for LRR and LASSO.
%
% Inputs:
%   - model
%     Either 'LRR' or 'LASSO'.
%
%   - result_dir
%     This refers to the directory in which the results were saved. Assumes a sub-directory
%     for each regression model (KRR, LRR, LASSO and RF).
%
%   - output_dir
%     This refers to the directory to save the collated results.
%
% Outputs:
%   - acc_<sample_size>_<model>
%     A mat file with a matrix of #folds/2 x #behaviours. Accuracy is in terms of
%     correlation averaged over the two split halves.
%
%   - icc_<sample_size>_<model>_Haufe
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances (Haufe) values from the two split halves.
%
%   - icc_<sample_size>_<model>_weights
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances (regression weights) values from the two split halves.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% setting up
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'predict_phenotypes','ChenOoi2023_ICCW', 'analysis', 'utilities');
addpath(genpath(project_code_dir));

num_behav = 1; % change to 39 for all behaviors
s_sizes = [800 400]; % change to [5260 3000 2000 800 400] for all samples
mkdirp(output_dir);

% load FC matrix
load(fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', ...
    'ChenTam2022_TRBPC', 'FC', 'FC_subjects_rs_all_score_mf.mat'));
FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);

% save accuracies, Haufe transform and original weights ICC in mat files
for s = 1:length(s_sizes)
    disp(['Now calculating for sample size : ', num2str(s_sizes(s))])
    N = num2str(s_sizes(s));
    % load fold information
    load(fullfile(result_dir, 'KRR', N, 'no_relative_5_fold_sub_list.mat'))
    icc_Haufe = zeros(126,num_behav);
    icc_weights = zeros(126,num_behav);
    acc = zeros(126,num_behav);
    
    % iterate over folds
    for i = 1:126
        y_pred_1 = [];
        y_pred_2 = [];
        % iterate over behaviors
        for j = 1:num_behav
            % get averaged accuracy and icc for regression weights
            r1 = load(fullfile(result_dir, model, N, strcat('behav_', num2str(j)), 'rng1', ...
                'params', strcat('fold_', num2str(i)), 'acc_test_all_score.mat'));
            acc1 = r1.acc_corr;
            beta1 = r1.beta;
            y_pred_1(:,end+1) = r1.y_pred_train;
            r2 = load(fullfile(result_dir, model, N, strcat('behav_', num2str(j)), 'rng1', ...
                'params', strcat('fold_', num2str(253-i)), 'acc_test_all_score.mat'));
            acc2 = r2.acc_corr;
            beta2 = r2.beta;
            y_pred_2(:,end+1) = r2.y_pred_train;
            
            % final accuracy
            acc(i,j) = (acc1+acc2)/2;
            % regression weight icc
            icc_weights(i,j) = CBIG_ICCW_ICC_1_1([beta1 beta2]);
        end
        
        % Haufe transform
        pfm1 = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(i).fold_index==0)',y_pred_1);
        pfm2 = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(253-i).fold_index==0)',y_pred_2);
        for j = 1:num_behav
            icc_Haufe(i,j) = CBIG_ICCW_ICC_1_1([pfm1(:,j) pfm2(:,j)]);
        end
    end
    
    % save results
    save(fullfile(output_dir, strcat('acc_', N, '_', model, '.mat')), 'acc');
    icc = icc_Haufe;
    save(fullfile(output_dir, strcat('icc_', N, '_', model, '_Haufe.mat')), 'icc');
    icc = icc_weights;
    save(fullfile(output_dir, strcat('icc_', N, '_', model, '_weights.mat')), 'icc');
end

rmpath(genpath(project_code_dir));
end