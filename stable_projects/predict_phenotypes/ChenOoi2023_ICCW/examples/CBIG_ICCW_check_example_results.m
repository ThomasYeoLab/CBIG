function CBIG_ICCW_check_example_results(out_dir)

% CBIG_ICCW_check_example_results(out_dir)
% This function checks if the generated example results are identical to
% the reference files.
%
% Input:
%   - out_dir
%     The output directory path saving results of example scripts
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% get directories
ref_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', ...
    'ChenOoi2023_ICCW', 'examples', 'ref_results');

% compare prediction accuracies of KRR
ref_acc = load(fullfile(ref_dir, 'final_result_2cog.mat'));
test_acc = load(fullfile(out_dir, 'final_result_2cog.mat'));

diff_acc = max(max(abs(ref_acc.optimal_acc - test_acc.optimal_acc)));
assert(diff_acc< 1e-5, sprintf('Difference in acc of KRR of: %f', diff_acc));

% compare the Haufe-transformed regression weights
ref_pfm = load(fullfile(ref_dir, 'pfm_icc.mat'));
test_pfm = load(fullfile(out_dir, 'pfm_icc.mat'));

diff_pfm1 = max(max(abs(ref_pfm.pfm1 - test_pfm.pfm1)));
assert(diff_acc< 1e-5, sprintf('Difference in pfm1 of: %f', diff_pfm1));
diff_pfm2 = max(max(abs(ref_pfm.pfm2 - test_pfm.pfm2)));
assert(diff_acc< 1e-5, sprintf('Difference in pfm2 of: %f', diff_pfm2));

% compare the icc of the Haufe-transformed regression weights
diff_icc = ref_pfm.icc - test_pfm.icc;
assert(diff_acc< 1e-5, sprintf('Difference in icc of: %f', diff_icc));

display('ChenOoi2023 example run successfully!')

end