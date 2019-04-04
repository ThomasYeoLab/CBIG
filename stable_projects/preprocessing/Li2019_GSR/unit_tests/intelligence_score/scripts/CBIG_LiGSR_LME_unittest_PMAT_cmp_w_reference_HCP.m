function CBIG_LiGSR_LME_unittest_PMAT_cmp_w_reference_HCP( all_stats_mat )

% CBIG_LiGSR_LME_unittest_PMAT_cmp_w_reference_HCP( all_stats_mat )
% 
% This function compare the variance component model results of fluid
% intelligence measure in HCP dataset that the user generated with the
% ground truth results.
% 
% Inputs:
%   - all_stats_mat
%     The full path of a .mat file. This file is the output of
%     CBIG_LiGSR_LME_unittest_PMAT_cmp2pipe_HCP.m
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ref_dir = fullfile('/mnt', 'eql', 'yeo1', 'CBIG_private_data', 'unit_tests', ...
    'stable_projects', 'preprocessing', 'Li2019_GSR');
ref_result = fullfile(ref_dir, 'intelligence_score', 'VarianceComponentModel', 'HCP', ...
    'ref_output', 'compare_2pipe', 'allstats_cmp2pipelines.mat');

ref = load(ref_result);
load(all_stats_mat);

assert(abs(ref.perc_improv - perc_improv)<1e-10, ...
    'Your percentage improvement differs from the reference.\n');
assert(abs(ref.m_jack - m_jack)<1e-10, ...
    'Your jackknife mean differs from the reference.\n');
assert(abs(ref.v_jack - v_jack)<1e-10, ...
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

