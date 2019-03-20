function [Params, lh_labels, rh_labels] = CBIG_MSHBM_example_wrapper(out_dir)

% [Params, lh_labels, rh_labels] = CBIG_MSHBM_example_wrapper(out_dir)
%
% This function will run examples sequentially
% 
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

%% Generating input example data
system(['${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Kong2019_MSHBM/examples/CBIG_MSHBM_create_example_input_data.sh ' out_dir]);


%% Example1: Group priors estimation
cd(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'step2_estimate_priors'));

project_dir = fullfile(out_dir, 'estimate_group_priors');
Params = CBIG_MSHBM_estimate_group_priors(project_dir,'fsaverage5','2','2','17','5');

%% Example2: Individual-level parcellations generation
cd(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'step3_generate_ind_parcellations'));

project_dir = fullfile(out_dir, 'generate_individual_parcellations');

% Generate individual parcellation for subject 1 using 2 sessions
[lh_labels, rh_labels] = CBIG_MSHBM_generate_individual_parcellation(project_dir,'fsaverage5','2','17','1','100','50');

