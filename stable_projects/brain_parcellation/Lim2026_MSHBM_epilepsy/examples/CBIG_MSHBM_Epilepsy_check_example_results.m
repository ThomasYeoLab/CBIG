function CBIG_MSHBM_Epilepsy_check_example_results(out_dir)

% ====================================
%% Help
% ====================================
%
% CBIG_MSHBM_Epilepsy_check_example_results
%
% Checks that the individual-specific parcellation and language laterality
% index (LI) produced by CBIG_MSHBM_Epilepsy_wrapper match the reference
% outputs stored in examples/results/. The reference outputs in examples/results/
% should not be overwritten; pass a separate out_dir to the wrapper.
%
% INPUTS:
%   out_dir: Path to the directory where the wrapper saved its outputs.
%
% OUTPUTS:
%   Prints a confirmation message if results match; raises an assertion
%   error with the maximum label difference if they do not.
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Load Reference and Test Parcellations
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
    'brain_parcellation', 'Lim2026_MSHBM_epilepsy', 'examples', 'results');

parcel_subpath = fullfile('ind_parcellation', 'test_set', ...
    'Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

ref  = load(fullfile(ref_dir, parcel_subpath));
test = load(fullfile(out_dir,  parcel_subpath));

% ====================================
%% Compare Labels
% ====================================

diff_lh = max(abs(double(test.lh_labels) - double(ref.lh_labels)));
diff_rh = max(abs(double(test.rh_labels) - double(ref.rh_labels)));

assert(diff_lh == 0, sprintf('LH labels differ from reference: max diff = %f', diff_lh));
assert(diff_rh == 0, sprintf('RH labels differ from reference: max diff = %f', diff_rh));

% ====================================
%% Compare Language Laterality Index
% ====================================

ref_LI  = load(fullfile(ref_dir, 'LI.mat'));
test_LI = load(fullfile(out_dir,  'LI.mat'));

diff_LI = abs(test_LI.LI_lang - ref_LI.LI_lang);

assert(diff_LI < 1e-6, sprintf('LI_lang differs from reference: diff = %f', diff_LI));

fprintf('Example results match reference.\n');

end
