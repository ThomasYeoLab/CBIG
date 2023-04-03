function CBIG_GradPar_check_example_results(out_dir)

% CBIG_GradPar_check_example_results(out_dir)
% This function checks if the generated example results are identical with
% the reference files.
%
% Input:
%   - out_dir: The output directory path saving results of example scripts (i.e. CBIG_GradPar_example_wrapper.m)
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

% check the results
ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes',...
'Kong2023_GradPar', 'examples', 'ref_results');

% compare prediction accuracies of KRR
for res = {'100', '200'}
    ref_acc = load(fullfile(ref_dir, 'KRR', res{1}, 'final_result_example_targets.mat'));
    test_acc = load(fullfile(out_dir, 'KRR', res{1}, 'final_result_example_targets.mat'));
    
    diff_acc = max(max(abs(ref_acc.optimal_acc - test_acc.optimal_acc)));
    assert(diff_acc< 1e-5, sprintf('maximum difference of KRR with resolution %s : %f', res{1}, diff_acc));
end
% compare prediction accuracies of LRR frac
for res = {'100', '200'}
    ref_acc = load(fullfile(ref_dir, 'LRR_frac', res{1}, 'results', 'optimal_acc', 'example_targets.mat'));
    test_acc = load(fullfile(out_dir, 'LRR_frac', res{1}, 'results', 'optimal_acc', 'example_targets.mat'));
    
    diff_acc = max(max(abs(ref_acc.acc_corr_test - test_acc.acc_corr_test)));
    assert(diff_acc< 1e-5, sprintf('maximum difference of LRR frac with resolution %s : %f', res{1}, diff_acc));
end

fprintf('Your example results are correct!\n')
