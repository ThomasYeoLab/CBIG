function p_value = CBIG_MMP_compute_HCP_permutation_p_value(ref_dir, behav_ind, outstem, perm_dir, seeds_total)

% p_value = CBIG_MMP_compute_HCP_permutation_p_value(ref_dir, behav_ind, outstem, perm_dir, seeds_total)
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
%   - seeds_total
%     Number of splits. Number of outerfolds and permutations are set to 10 and 10000 respectively.
%
% Outputs:
%   - p-value
%     A vector. P-value for behaviour in behav_ind.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add utilities path
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'HCP', 'utilities'))

%% compare with permutation: MODIFY HERE IF YOU CHANGED FOLDS / PERMUTATIONS
perm_per_file = 10000;
N_folds = 10;

for i = 1:length(behav_ind)
    N_worse = 0;
    behav = behav_ind(i);
    
    % load reference
    acc_vec = CBIG_MMP_HCP_read_model_results(outstem, ref_dir, seeds_total, N_folds, behav_ind, 'corr', 0);
    for curr_seed = 1:seeds_total
        seed_str = strcat('seed_', num2str(curr_seed));
        seed_alloc = [ ((curr_seed-1)*N_folds + 1):curr_seed*N_folds ];
        ref = mean(acc_vec(i,seed_alloc));
        % load permutations
        load(fullfile(perm_dir, seed_str, ['/acc_score' num2str(behav) '_allFolds_permStart1.mat']));
        curr_acc = squeeze(mean(stats_perm.corr,1));
        curr_N_worse = sum(curr_acc > ref);
        N_worse = N_worse + curr_N_worse;
    end
    p_value(i) = (1+N_worse)/(1+seeds_total*perm_per_file);
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'HCP', 'utilities'))
