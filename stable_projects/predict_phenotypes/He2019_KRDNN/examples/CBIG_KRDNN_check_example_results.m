function  CBIG_KRDNN_check_example_results(out_dir, opt_flag)

% CBIG_KRDNN_check_example_results(out_dir)
% This function checks if the generated example results are identical with
% the reference files.
%
% Input:
%   - out_dir: The output directory path saving results of example scripts (i.e. CBIG_MSHBM_example_wrapper.m)
%   - opt_flag: bool, optional (default is true), flag indicate checking CV or TVT results, true is CV, false 
%               is TVT
%
% Written by He Tong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin > 1
  flag_cv = opt_flag;
else
  flag_cv = true;
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

% check the results
if flag_cv
    ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
        'predict_phenotypes', 'He2019_KRDNN', 'examples', ...
        'results', 'KR_CV_matlab_r2018b');
    
    ref_output = load(fullfile(ref_dir, 'final_result.mat'));
    test_output = load(fullfile(out_dir, 'final_result.mat'));

    % check accuracy
    diff_acc = sum(sum(abs(ref_output.optimal_acc - test_output.optimal_acc)));
    assert(diff_acc < 1e-6, sprintf('optimal_acc difference: %f', diff_acc));
else
    ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
        'predict_phenotypes', 'He2019_KRDNN', 'examples', ...
        'results', 'KR_TVT');
    
    ref_output = load(fullfile(ref_dir, 'final_result.mat'));
    test_output = load(fullfile(out_dir, 'final_result.mat'));

    % check correlation
    diff_corr = sum(abs(ref_output.optimal_corr - test_output.optimal_corr));
    assert(diff_corr < 1e-6, sprintf('optimal_corr difference: %f', diff_corr));

    % check MAE of age prediction
    diff_mae = abs(ref_output.age_mae - test_output.age_mae);
    assert(diff_mae < 1e-6, sprintf('age_mae difference: %f', diff_mae));

    % check accuracy of sex prediction
    diff_acc = abs(ref_output.sex_accuracy - test_output.sex_accuracy);
    assert(diff_acc < 1e-6, sprintf('sex_accuracy difference: %f', diff_acc));
end

fprintf('Your example results are correct!\n')