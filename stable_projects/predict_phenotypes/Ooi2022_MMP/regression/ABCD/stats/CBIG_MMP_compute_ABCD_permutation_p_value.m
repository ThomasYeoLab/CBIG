function p_value = CBIG_MMP_compute_ABCD_permutation_p_value(ref_dir, behav_ind, outstem, perm_dir, perm_per_file)

% p_value = CBIG_MMP_compute_ABCD_permutation_p_value(ref_dir, behav_ind, outstem, perm_dir, perm_per_file)
%
% Adapted from CBIG_TRBPC_compute_permutation_p_value.
% This function computes of p-value of permutation test.
%
% Inputs:
%   - ref_dir
%     Directory of the regression results.
%
%   - behav_ind
%     Behaviourial indices to calculate p-value file.
%
%   - outstem
%     Name of regression results folder (e.g. KRR_features_ct).
%
%   - perm_dir
%     Directory of permutation results.
%
%   - perm_per_file
%     Number of permutations in permutation test.
%
% Outputs:
%   - p-value
%     A vector. P-value for behaviour in behav_ind.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add utilities path
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'))

%% compare with permutation
for i = 1:length(behav_ind)
    N_worse = 0;
    behav = behav_ind(i);
    
    % load reference
    acc_vec = CBIG_MMP_ABCD_read_model_results(outstem, ref_dir, 120, behav_ind, 'corr', 0);
    ref = mean(acc_vec(i,:));
    % load permutations
    load(fullfile(perm_dir, ['acc_score' num2str(behav) '_allFolds_permStart1.mat']));
    curr_acc = squeeze(mean(stats_perm.corr,1));
    curr_N_worse = sum(curr_acc > ref);
    N_worse = N_worse + curr_N_worse;
    
    p_value(i) = (1+N_worse)/(1+perm_per_file);
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'utilities'))
