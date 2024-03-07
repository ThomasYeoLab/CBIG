function [Params, lh_labels, rh_labels] = CBIG_MSHBM_example_wrapper(out_dir)

% [Params, lh_labels, rh_labels] = CBIG_MSHBM_example_wrapper(out_dir)
%
% This function will run examples sequentially
% 
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

%% Example 101: A simplified way to generate individual parcellation for a single subject
CBIG_MSHBM_example_single_subject(fullfile(out_dir, 'single_subject_parcellation'));


%% Generating input example data
cmd = '${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Kong2019_MSHBM';
cmd = [cmd '/examples/CBIG_MSHBM_create_example_input_data.sh ' out_dir];
system(cmd);


%% Example1: Group priors estimation
cd(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'step2_estimate_priors'));

project_dir = fullfile(out_dir, 'estimate_group_priors');
Params = CBIG_MSHBM_estimate_group_priors(project_dir, 'fsaverage5', '2', '2', '17', 'max_iter', '5');

%% Example2: Individual-level parcellations generation
cd(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', ...
    'step3_generate_ind_parcellations'));

project_dir = fullfile(out_dir, 'generate_individual_parcellations');

% Generate individual parcellation for subject 1 using 2 sessions
[lh_labels, rh_labels] = CBIG_MSHBM_generate_individual_parcellation(project_dir,'fsaverage5','2','17','1','100','50');