function CBIG_LiGSR_LRR_replication_cmp_w_reference_GSP( final_result )

% CBIG_LiGSR_LRR_replication_cmp_w_reference_GSP( final_result )
% 
% This function compares the unit test results generated using
% CBIG_LiGSR_LRR_replication_cmp2pipe_GSP against the ground truth results
% to determine if all the related linear ridge regression function are
% working accurately.
% 
% Inputs:
%   - final_result
%     The full path of a .mat file. This file is the output of
%     CBIG_LiGSR_LRR_replication_cmp2pipe_GSP.m
%     It should contain three 20x24 matrices: 
%     acc_GSR_mean: the mean accuracy of each random split of 24 behavioral
%                   and demographic measures with global signal regression.
%     acc_Baseline_mean: the mean accuracy of each random split of 24
%                        behavrioal and demographic measures with baseline
%                        fMRI preprocessing.
%     mean_acc_dif: the difference between "acc_GSR_mean" and
%                   "acc_Baseline_mean".
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ref_dir = fullfile(getenv('CBIG_REPDATA_DIR'), ...
    'stable_projects', 'preprocessing', 'Li2019_GSR');
ref_result = fullfile(ref_dir, 'LinearRidgeRegression', ...
    'GSP', 'ref_output', 'compare_2pipe', 'final_result.mat');

ref = load(ref_result);
load(final_result);

assert(max(max(abs(ref.acc_GSR_mean - acc_GSR_mean)))<1e-15, ...
    'Your accuracies with GSR differ from the reference.');
assert(max(max(abs(ref.acc_Baseline_mean - acc_Baseline_mean)))<1e-15, ...
    'Your baseline accuracies differ from the reference.');
assert(max(max(abs(ref.mean_acc_dif - mean_acc_dif)))<1e-15, ...
    'Your accuracy differences between pipelines deviate from the reference.\n')

fprintf('Your results replicated the reference results.\n')

end

