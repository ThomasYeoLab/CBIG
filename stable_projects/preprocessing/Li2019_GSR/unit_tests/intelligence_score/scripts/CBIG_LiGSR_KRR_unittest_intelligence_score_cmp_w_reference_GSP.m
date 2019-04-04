function CBIG_LiGSR_KRR_unittest_intelligence_score_cmp_w_reference_GSP( final_result )

% CBIG_LiGSR_KRR_unittest_intelligence_score_cmp_w_reference_GSP( final_result )
% 
% This function compares the accuracy of predicting Shipley Vocabulary and
% WAIS - Matrix Reasoning measures in GSP dataset using kernel ridge
% regression against the ground truth results.
% 
% Inputs:
%   - final_result
%     The full path of a .mat file. This file is the output of
%     CBIG_LiGSR_KRR_unittest_intelligence_score_cmp2pipe_GSP.m
%     It should contain three 20x2 matrices: 
%     acc_GSR_mean: the mean accuracy of each random split of the 2
%                   behavioral measures with global signal regression.
%     acc_Baseline_mean: the mean accuracy of each random split of the 2
%                        behavrioal measures with baseline fMRI
%                        preprocessing.
%     mean_acc_dif: the difference between "acc_GSR_mean" and
%                   "acc_Baseline_mean".
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ref_dir = fullfile('/mnt', 'eql', 'yeo1', 'CBIG_private_data', 'unit_tests', ...
    'stable_projects', 'preprocessing', 'Li2019_GSR');
ref_result = fullfile(ref_dir, 'intelligence_score', 'KernelRidgeRegression', ...
    'GSP', 'ref_output', 'compare_2pipe', 'final_result.mat');

ref = load(ref_result);
load(final_result);

assert(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean)))<1e-10, ...
    'Your accuracies with GSR differ from the reference.');
assert(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean)))<1e-10, ...
    'Your baseline accuracies differ from the reference.');
assert(max(max(abs(ref.mean_acc_dif - mean_acc_dif)))<1e-10, ...
    'Your accuracy differences between pipelines deviate from the reference.')

fprintf('Your results replicated the reference results.\n')

end

