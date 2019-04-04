function CBIG_LiGSR_LME_replication_cmp_w_reference_HCP( all_stats_mat )

% CBIG_LiGSR_LME_replication_cmp_w_reference_HCP( all_stats_mat )
% 
% This function compares the variance component model results of HCP
% dataset that the user generated to replicate Li et al., 2019 with the
% ground truth results.
% 
% Inputs:
%   - all_stats_mat
%     The full path of a .mat file. This file is the output of
%     CBIG_LiGSR_LME_replication_cmp2pipe_HCP.m
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ref_dir = fullfile('/mnt', 'eql', 'yeo1', 'CBIG_private_data', 'replication', ...
    'stable_projects', 'preprocessing', 'Li2019_GSR');
ref_result = fullfile(ref_dir, 'VarianceComponentModel', 'HCP', ...
    'ref_output', 'compare_2pipe', 'allstats_cmp2pipelines.mat');

ref = load(ref_result);
load(all_stats_mat);

assert(abs(ref.perc_improv - perc_improv)<1e-15, ...
    'Your percentage improvement differs from the reference.\n');
assert(abs(ref.m_jack - m_jack)<1e-15, ...
    'Your jackknife mean differs from the reference.\n');
assert(abs(ref.v_jack - v_jack)<1e-15, ...
    'Your jackknife variance differs from the reference.\n')
assert(abs(ref.IQR_pos - IQR_pos)<1e-15, ...
    'Your "#traits whose IQR > 0" differs from the reference.\n')
assert(abs(ref.IQR_neg - IQR_neg)<1e-15, ...
    'Your "#traits whose IQR < 0" differs from the reference.\n')
assert(abs(ref.med_pos - med_pos)<1e-15, ...
    'Your "#traits whose median > 0" differs from the reference.\n')
assert(abs(ref.med_pos - med_pos)<1e-15, ...
    'Your "#traits whose median < 0" differs from the reference.\n')

fprintf('Your results replicated the reference results.\n')

end

