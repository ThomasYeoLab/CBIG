function  CBIG_MM_check_example_results(out_dir)

% CBIG_MM_check_example_results(out_dir)
% This function checks if the generated example results are identical with
% the reference files.
%
% Input:
%   - out_dir: The output directory path saving results of example scripts (i.e. CBIG_MSHBM_example_wrapper.m)
%
% Written by He Tong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

% check the results
ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
    'predict_phenotypes', 'He2022_MM', 'examples', ...
    'results');

ref_output = load(fullfile(ref_dir, 'final_result.mat'));
test_output = load(fullfile(out_dir, 'final_result.mat'));

% check accuracy
diff_acc = sum(sum(abs(ref_output.optimal_acc - test_output.optimal_acc)));
assert(diff_acc < 1e-6, sprintf('optimal_acc difference: %f', diff_acc));


fprintf('Your example results are correct!\n')