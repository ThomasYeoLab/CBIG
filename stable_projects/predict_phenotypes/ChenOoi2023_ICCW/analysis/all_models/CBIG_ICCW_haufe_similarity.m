function CBIG_ICCW_haufe_similarity(result_dir, output_dir)

% function CBIG_ICCW_haufe_similarity(result_dir, output_dir)
%
% This function computes the similarity between the Haufe transformed weights across
% regression models.
%
% Inputs:
%   - result_dir
%     This refers to the directory in which the results were saved. Assumes a sub-directory
%     for each regression model (KRR, LRR, LASSO and RF).
%
%   - output_dir
%     This refers to the directory to save the collated results.
%
% Outputs:
%   - sim_behav
%     A 4x4 matrix containing the similarity values between KRR, LRR, LASSO and RF.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% setting up
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'ChenOoi2023_ICCW', 'analysis', 'utilities');
addpath(genpath(project_code_dir));

sim_behav = zeros(4,4,36);
% load FC matrix
load(fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', ...
    'ChenTam2022_TRBPC', 'FC', 'FC_subjects_rs_all_score_mf.mat'));
% normalize FC
FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);

% load y and folds
load(fullfile(result_dir, 'KRR', num2str(5260), 'final_result_all_score.mat'), 'y_pred_train');
load(fullfile(result_dir, 'KRR', num2str(5260), 'no_relative_5_fold_sub_list.mat'));
y_pred_KRR = y_pred_train;

% iterate over folds
for i = 1:252
    y_pred_LRR = [];
    y_pred_LASSO = [];
    
    % iterate over behaviors
    for j = 1:36
        % extract LRR and LASSO y pred values to compute Haufe
        load(fullfile(result_dir, 'LRR', num2str(5260), strcat(behav_' num2str(j)), 'rng1', 'params', ...
            strcat('fold_', num2str(i)), 'acc_test_all_score.mat'));
        y_pred_LRR(:,end+1) = y_pred_train;
        
        load(fullfile(result_dir, 'LASSO', num2str(5260), strcat(behav_' num2str(j)), 'rng1', 'params', ...
            strcat('fold_', num2str(i)), 'acc_test_all_score.mat'));
        y_pred_LASSO(:,end+1) = y_pred_train;
        
        % read Haufe values from csv
        fi_csv = csvread(fullfile(result_dir, 'RF', num2str(5260), strcat(behav_' num2str(j)), ...
            'rng0', strcat('fold_', num2str(i)), strcat('ABCD_RF_behav_', num2str(j), '_fi.csv')), 1,0);
        pfm_RF(:,j) = fi_csv(:,4);
    end
    % recompute Haufe values
    pfm_KRR = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(i).fold_index==0)',y_pred_KRR{i});
    pfm_LRR = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(i).fold_index==0)',y_pred_LRR);
    pfm_LASSO = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(i).fold_index==0)',y_pred_LASSO);
    
    %  compare to KRR
    sim_behav(2,1,:) = diag(corr(pfm_KRR, pfm_LRR));
    sim_behav(3,1,:) = diag(corr(pfm_KRR, pfm_LASSO));
    sim_behav(4,1,:) = diag(corr(pfm_KRR, pfm_RF));
    % compare to LRR
    sim_behav(3,2,:) = diag(corr(pfm_LRR, pfm_LASSO));
    sim_behav(4,2,:) = diag(corr(pfm_LRR, pfm_RF));
    % compare to LASSO
    sim_behav(4,3,:) = diag(corr(pfm_LASSO, pfm_RF));
    
end
% save final results
save(fullfile(output_dir, 'sim_behav_Haufe.mat'), 'sim_behav')

rmpath(genpath(project_code_dir));
end