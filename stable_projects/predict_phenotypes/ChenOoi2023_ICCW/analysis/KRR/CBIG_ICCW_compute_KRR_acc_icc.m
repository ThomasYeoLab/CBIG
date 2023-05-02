function CBIG_ICCW_compute_KRR_acc_icc(krr_result_dir, krr_weights_folder, output_dir)

% function CBIG_ICCW_compute_KRR_acc_icc(krr_result_dir, krr_weights_folder, output_dir)
%
% This function collates the accuracy results, and computes the ICC for the Haufe
% transformed weights, as well as the original regression weights for KRR
%
% Inputs:
%   - krr_result_dir
%     This refers to the directory in which the KRR results were saved.
%
%   - krr_weights_folder
%     This refers to the sub-directory within krr_results_dir in which the
%     KRR weights were saved.
%
%   - output_dir
%     This refers to the directory to save the collated results.
%
% Outputs:
%   - acc_<sample_size>_KRR
%     A mat file with a matrix of #folds/2 x #behaviours. Accuracy is in terms of
%     correlation averaged over the two split halves.
%
%   - icc_<sample_size>_KRR_Haufe
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances (Haufe) values from the two split halves.
%
%   - icc_<sample_size>_KRR_weights
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances (regression weights) values from the two split halves.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% setting up
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'predict_phenotypes', 'ChenOoi2023_ICCW', 'analysis', 'utilities');
addpath(genpath(project_code_dir));

num_behav = 1; % change to 39 for all behaviors
s_sizes = [800 400]; % change to [5260 3000 2000 800 400] for all samples
mkdirp(output_dir)

% load FC matrix
load(fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', ...
    'ChenTam2022_TRBPC', 'FC', 'FC_subjects_rs_all_score_mf.mat'));
% normalize FC
FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);

% save accuracies in mat files
for samples = 1:length(s_sizes)
    disp(['Now calculating for sample size : ', num2str(s_sizes(samples))])
    load(fullfile(krr_result_dir,num2str(s_sizes(samples)), 'final_result_all_score.mat'))
    load(fullfile(krr_result_dir,num2str(s_sizes(samples)), 'no_relative_5_fold_sub_list.mat'))
    % save accuracies in mat file
    acc = zeros(126,39);
    for i = 1:126
        acc(i,:) = (optimal_acc(i,:) + optimal_acc(253-i))/2;
    end
    save(fullfile(output_dir, strcat('acc_', num2str(s_sizes(samples)), '_KRR.mat')), 'acc');
    
    % save Haufe transform ICCs in mat files
    icc = zeros(126,num_behav);
    for i = 1:126
        pfm1 = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(i).fold_index==0)',y_pred_train{i});
        pfm2 = CBIG_ICCW_cov_matrix(FC_all(:,sub_fold(253-i).fold_index==0)',y_pred_train{253-i});
        for j = 1:num_behav
            icc(i,j) = CBIG_ICCW_ICC_1_1([pfm1(:,j) pfm2(:,j)]);
        end
    end
    save(fullfile(output_dir, strcat('icc_', num2str(s_sizes(samples)), '_KRR_Haufe.mat')), 'icc');
    
    
    % save regression weight ICCs in mat files
    icc = zeros(126,num_behav);
    for j = 1:num_behav
        load(fullfile(krr_result_dir, num2str(s_sizes(samples)), krr_weights_folder, ...
            strcat('weights_score', num2str(j), '_all_folds.mat')));
        for i = 1:126
            icc(i,j) = CBIG_ICCW_ICC_1_1([weights_all_folds(:,i), weights_all_folds(:,253-i)]);
        end
    end
    save(fullfile(output_dir, strcat('icc_', num2str(s_sizes(samples)), '_KRR_weights.mat')), 'icc')
    
end

rmpath(genpath(project_code_dir));
end
