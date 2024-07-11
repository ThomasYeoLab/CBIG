function  CBIG_MMM_check_example_results(ref_output_dir, output_classical_dir)

% CBIG_MMM_check_example_results(out_dir)
% This function checks if the generated example results are identical with
% the reference files.
%
% Input:
%   - ref_output_dir: The output directory path saving refernce results
%   - output_classical_dir: The output directory path saving of KRR classical
%
% Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% get reference output folder
ref_output = load(fullfile(ref_output_dir, 'output_KRR_classical_exp_test', 'final_result', ...
    'krr_classical_res.mat'));
test_output = load(fullfile(output_classical_dir, 'final_result', 'krr_classical_res.mat'));

% check KRR classical accuracy
diff_cod = sum(sum(sum(abs(ref_output.meta_cod - test_output.meta_cod))));
assert(diff_cod < 1e-6, sprintf('KRR classical meta_cod difference: %f', diff_cod));
diff_cor = sum(sum(sum(abs(ref_output.meta_cor - test_output.meta_cor))));
assert(diff_cor < 1e-6, sprintf('KRR classical meta_cor difference: %f', diff_cor));

fprintf('Your example results are correct!\n')