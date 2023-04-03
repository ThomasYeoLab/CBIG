function CBIG_hMRF_example_wrapper(output_dir)
% CBIG_hMRF_example_wrapper(output_dir)
%
% This is the wrapper function to set up and run the clustering algorithm based on the input arguments.
% 
% Input
%   - output_dir: (string) 
%     ABSOLUTE path to the directory to which the output results will be saved.
%
% Example
%   - CBIG_hMRF_example_wrapper('/home/user/example_results/output_folder')
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
hMRF_code_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Yan2023_homotopic', 'code');
hMRF_example_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation',...
    'Yan2023_homotopic', 'examples');
addpath(genpath(hMRF_code_dir));

% remove the output folder if exist
if(exist(output_dir, 'dir'))
    disp('Output folder already exists. Removing output folder...');
    rmdir(output_dir, 's')
end
mkdir(output_dir);

%% step1: create input subject list
lh_fmri_fullpath_txt = fullfile(output_dir, 'lh_GSP_subject_fullpath.csv');
rh_fmri_fullpath_txt = fullfile(output_dir, 'rh_GSP_subject_fullpath.csv');
censor_file_fullpath = fullfile(output_dir, 'censor_full_path.csv');
cmd = ['sh ' fullfile(hMRF_example_dir, 'CBIG_hMRF_create_2subject_fullpaths.sh') ' ' lh_fmri_fullpath_txt...
    ' ' rh_fmri_fullpath_txt ' ' censor_file_fullpath];
system(cmd); pause(5);

%% step2: create premultiplied matrix
CBIG_hMRF_generate_premultiplied_matrix(output_dir, lh_fmri_fullpath_txt, rh_fmri_fullpath_txt, 'fsaverage6',...
censor_file_fullpath);

%% step3: run parcellation algorithm
parc_outdir = fullfile(output_dir, 'parcellation');
mkdir(parc_outdir);
cd(fullfile(hMRF_code_dir, 'step2_generate_parcellation'));

results = CBIG_hMRF_wrapper_generate_homotopic_parcellation(fullfile(output_dir, 'premultiplied_matrix_single.mat'),...
  parc_outdir,'start_seed', 835, 'num_rand_inits', 1, 'num_cluster', 100, 'num_iterations', 3, 'initial_c', 100,...
   'initial_d', 10, 'k', 15, 'w_xyz', 1500, 'mesh_type', 'fsaverage6', 'decrease_c', true, 'decrease_d', true,...
   'increase_tau', true, 'premature_stopping', false);

% visualize final parcellation output
% your figure should look very similar to the figure under `./ref_results/parcellation_seed_835`
CBIG_DrawSurfaceMaps(results.lh_label, results.rh_label, 'fsaverage6', 'inflated', 0, 100, colorcube);
rmpath(genpath(hMRF_code_dir));
end
