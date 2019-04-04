function CBIG_LiGSR_KRR_replication_cmp_w_reference_HCP( final_result )

% CBIG_LiGSR_KRR_replication_cmp_w_reference_HCP( final_result )
% 
% This function compares the unit test results generated using
% CBIG_LiGSR_KRR_replication_cmp2pipe_HCP against the ground truth results
% to determine if all the related kernel ridge regression function are
% working accurately.
% 
% Inputs:
%   - final_result
%     The full path of a .mat file. This file is the output of
%     CBIG_LiGSR_KRR_replication_cmp2pipe_HCP.m
%     It should contain three 20x60 matrices: 
%     acc_GSR_mean: the mean accuracy of each random split of 60 behavioral
%                   and demographic measures with global signal regression.
%     acc_Baseline_mean: the mean accuracy of each random split of 60
%                        behavrioal and demographic measures with baseline
%                        fMRI preprocessing.
%     mean_acc_dif: the difference between "acc_GSR_mean" and
%                   "acc_Baseline_mean".
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ref_dir = fullfile('/mnt', 'eql', 'yeo1', 'CBIG_private_data', 'replication', ...
    'stable_projects', 'preprocessing', 'Li2019_GSR');
ref_result = fullfile(ref_dir, 'KernelRidgeRegression', ...
    'HCP', 'ref_output', 'compare_2pipe', 'final_result.mat');

ref = load(ref_result);
load(final_result);

assert(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean)))<1e-15, ...
    'Your accuracies with GSR differ from the reference.');
assert(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean)))<1e-15, ...
    'Your baseline accuracies differ from the reference.');
assert(max(max(abs(ref.mean_acc_dif - mean_acc_dif)))<1e-15, ...
    'Your accuracy differences between pipelines deviate from the reference.')

fprintf('Your results replicated the reference results.\n')

end

