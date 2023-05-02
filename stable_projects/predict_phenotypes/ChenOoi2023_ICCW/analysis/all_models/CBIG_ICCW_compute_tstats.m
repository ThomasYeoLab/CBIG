function CBIG_ICCW_compute_tstats(krr_result_dir, output_dir)

% function CBIG_ICCW_compute_KRR_acc_icc(krr_result_dir, output_dir)
%
% This function computes original regression weights for KRR. The original weights
% are a #features x #folds matrix which is saved into a mat file. Please ensure
% that KRR has completed running for all sample sizes first.
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
%   - icc_<sample_size>KRR_weights
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances values from the two split halves.
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
features_norm = FC_all';

% iterate over each sample size
for s = 1:length(s_sizes)
    disp(['Processing for sample size ', num2str(s_sizes(s))])
    curr_s = num2str(s_sizes(s));
    
    % load y and folds
    load(fullfile(krr_result_dir, curr_s, 'y_all_score.mat'));
    load(fullfile(krr_result_dir, curr_s, 'no_relative_5_fold_sub_list.mat'));
    
    % compute ICC of tstats for each fold
    icc = zeros(126,num_behav);
    for i = 1:126
        ind1 = sub_fold(i).fold_index == 1;
        ind2 = sub_fold(253-i).fold_index == 1;
        
        f1 = features_norm(ind1,:);
        f2 = features_norm(ind2,:);
        
        load(fullfile(krr_result_dir, curr_s, 'y', strcat('fold_', num2str(i)), 'y_regress_all_score.mat'));
        y1 = y_resid(ind1,:);
        load(fullfile(krr_result_dir, curr_s, 'y', strcat('fold_', num2str(253-i)), 'y_regress_all_score.mat'));
        y2 = y_resid(ind2,:);
        
        % compute t stats
        pfm1 = corr(f1,y1);
        pfm2 = corr(f2,y2);
        t1 = pfm1./(1-pfm1.^2);
        t2 = pfm2./(1-pfm2.^2);
        for j = 1:num_behav
            icc(i,j) = CBIG_ICCW_ICC_1_1([t1(:,j), t2(:,j)]);
        end
    end
    save(fullfile(output_dir, strcat('icc_', curr_s, '_tstats.mat')), 'icc')
end

rmpath(genpath(project_code_dir));
end
