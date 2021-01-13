function  CBIG_MSHBM_check_example_results(out_dir)

% CBIG_MSHBM_check_example_results(out_dir)
% This function checks if the generated example results are identical with
% the reference files.
%
% Input:
%   - out_dir: The output directory path saving results of example scripts (i.e. CBIG_MSHBM_example_wrapper.m)
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

% check the results
ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
    'brain_parcellation', 'Kong2019_MSHBM', 'examples', 'results');

% compare estimated group priors
ref_Params = load(fullfile(ref_dir, 'estimate_group_priors', 'priors', 'Params_Final.mat'));
test_Params = load(fullfile(out_dir, 'estimate_group_priors', 'priors', 'Params_Final.mat'));

diff_mu = max(max(abs(ref_Params.Params.mu - test_Params.Params.mu)));
assert(diff_mu < 1e-6, sprintf('maximum Params.mu difference: %f',diff_mu));

% compare estimated individual parcellation
ref_labels = load(fullfile(ref_dir, 'generate_individual_parcellations',...
    'ind_parcellation', 'Ind_parcellation_MSHBM_sub1_w100_MRF50.mat'));
test_labels = load(fullfile(out_dir, 'generate_individual_parcellations',...
    'ind_parcellation', 'test_set','Ind_parcellation_MSHBM_sub1_w100_MRF50.mat'));
diff_labels = max(abs([test_labels.lh_labels; test_labels.rh_labels] - [ref_labels.lh_labels; ref_labels.rh_labels]));
assert(diff_labels == 0, sprintf('maximum labels difference: %f',diff_labels));

fprintf('Your example results are correct!\n')
