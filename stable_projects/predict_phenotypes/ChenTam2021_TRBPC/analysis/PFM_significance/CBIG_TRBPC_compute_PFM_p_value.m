function p_value = CBIG_TRBPC_compute_PFM_p_value(PFM_network, perm_dir, behav_cat, perm_per_file, perm_total)

% p_value = CBIG_TRBPC_compute_PFM_p_value(ref_PFM, perm_dir, behav_cat, perm_per_file, perm_total)
%
% This function computes of p-value of predictive-feature matrix (PFM)
% edges by comparing the reference PFM with PFM in permutation test
%
% Inputs:
%   - PFM_network
%     the reference PFM averaged in each network block
%
%   - perm_dir
%     Directory of the PFM permutation results
%
%   - behav_cat
%     Path of the .mat behavior categories file. Inside the file a varible
%     named "behav_cat" is stored. The behav_cat variable is a cell array
%     of length # behavior categories. Inside each cell is a vector
%     indicate the indices for all behavior categories. For example, if
%     behav_cat{1} = [1,3,5] and  behav_cat{2} = [2,4,6], it means
%     behaviors 1,3,5 are in category 1 (e.g. cognition) and
%     behaviors 2,4,6 are in category 2 (e.g. mental health)
%
%   - perm_per_file
%     Number of permutations for each file in perm_dir
%
%   - perm_total
%     Total number of permutations
%
% Outputs:
%   - p-value
%     A vector. P-value for each edge in the ref_PFM
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% compute PFM avearged within each behavior category
N_edge = size(PFM_network,1);
PFM_mean_cat = zeros(N_edge,length(behav_cat));
N_cat = length(behav_cat);
for i = 1:N_cat
    PFM_mean_cat(:,i) = mean(PFM_network(:,behav_cat{i}),2);
end

%% compare with permutation
perm_seed = 1:perm_per_file:perm_total;
N_worse = zeros(N_edge,N_cat);
for i = 1:length(perm_seed)
    load(fullfile(perm_dir, ['/PFM_score_perm_start' num2str(perm_seed(i)) '.mat']),'PFM_score_perm');
    for curr_cat = 1:N_cat
        PFM_curr_cat = squeeze(mean(PFM_score_perm(:,behav_cat{curr_cat},:),2));
        curr_N_worse = sum(abs(PFM_curr_cat) > abs(PFM_mean_cat(:,curr_cat)),2);
        N_worse(:,curr_cat) = N_worse(:,curr_cat) + curr_N_worse;
    end
end

p_value = (1+N_worse)/(1+length(perm_seed)*perm_per_file);
end