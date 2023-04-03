function CBIG_hMRF_generate_premultiplied_matrix(output_dir, lh_fmri_fullpath_txt, rh_fmri_fullpath_txt, mesh_type,...
    censor_file_fullpath)

% CBIG_hMRF_generate_premultiplied_matrix(output_dir, lh_fmri_fullpath_txt, rh_fmri_fullpath_txt, mesh_type,...
% censor_file_fullpath)
% This function generates the premultiplied matrix with the given list of fmri input data.
%
% Input:
%   - output_dir: (char string)
%     The path to which the output is written. If directory does NOT exist, it will be created.
%
%   - <lh/rh>_fmri_fullpath_txt: (char string)
%     The path to the txt file that contains the list of fmri files for each subject, of <lh/rh> respectively.
%     Accepted fMRI file format: 'nii.gz'.
%     Each line should contain all fmri files of the current scan, separated by a single space.
% 
%   - mesh_type: (char string)
%     By default, 'fsaverage6'.
%
% Optional input:
%
%   - censor_file_fullpath: (char string)
%     The path to the txt file that contains the list of censor files for each subject, of right hemisphere.
%     Each censor file is a binary txt file, of length of total number of frames per run per subject.
%     Each line in this file should contain all censor files of the current scan, separated by a single space.
%
% Output:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate output files (can be REMOVED once the final_PMM matrix had been generated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - ${output_dir}/time_mat/<lh/rh>_time_mat.mat (matrix)
%     The concatenated and normalized fMRI time courses of NxT dimension. N is the number
%     of vertices of the given input surface mesh. T is the total number of time points in the resultant 
%     concatenated time courses.   
%   - ${output_dir}/mult_mat/<lh/rh>_multi_mat.mat (matrix)
%     The NxN premultipled matrix for each hemisphere, respectively. N is the number of vertices on the given surface
%     mesh, excluding the medial wall vertices.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - final_PMM (matrix)
%     This matrix is saved under fullfile(output_dir, 'premultiplied_matrix_single.mat').
%     It's the final premultiplied matrix across both hemispheres, with the left and right hemisphere matrices
%     occupying diagonal positions while all other entries are zero.
%
% Example:
%   CBIG_hMRF_generate_premultiplied_matrix('./output', './input/subjects_fullpath.txt', 'fsaverage6',...
%   './censor_files_fullpath.txt');
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~exist(output_dir, 'dir'))
    mkdir(output_dir);
end

if(~exist('censor_file_fullpath', 'var'))
    censor_file_fullpath = '';
end

mkdir(fullfile(output_dir,'time_data'));
mkdir(fullfile(output_dir,'mult_mat'));

lh_time_mat_file = fullfile(output_dir,'time_data','lh_time_matrix.mat');
rh_time_mat_file = fullfile(output_dir,'time_data','rh_time_matrix.mat');

% skip if input data already exists
if(exist(lh_time_mat_file, 'file') & exist(rh_time_mat_file, 'file'))
    disp('Output concatenated fMRI time matrices for both hemispheres already exists.');
else
    CBIG_hMRF_build_time_matrix(lh_fmri_fullpath_txt, rh_fmri_fullpath_txt,...
        mesh_type, lh_time_mat_file, rh_time_mat_file, censor_file_fullpath);
end

%% concatenate and process the premultiplied matrix if non-existent
lh_mult_file = fullfile(output_dir, 'mult_mat', 'lh_mult_matrix.mat');
rh_mult_file = fullfile(output_dir, 'mult_mat', 'rh_mult_matrix.mat');

if(exist(lh_mult_file, 'file') & exist(rh_mult_file, 'file'))
    disp('Output premultiplied fMRI time matrices for both hemispheres already exists.');
    load(lh_mult_file, 'lh_cov_mat');
    load(rh_mult_file, 'rh_cov_mat', 'dim');
else
    [lh_cov_mat, rh_cov_mat, dim] = CBIG_hMRF_build_prod_matrix(lh_time_mat_file, rh_time_mat_file,...
        lh_mult_file, rh_mult_file);
end

lh_mesh = CBIG_ReadNCAvgMesh('lh', mesh_type, 'inflated', 'cortex');
rh_mesh = CBIG_ReadNCAvgMesh('rh', mesh_type, 'inflated', 'cortex');
lh_mask = (lh_mesh.MARS_label == 2); % 1 for medial wall, 2 for cortex
rh_mask = (rh_mesh.MARS_label == 2); 

lh_full_cov_mat(lh_mask, lh_mask) = single(lh_cov_mat);
rh_full_cov_mat(rh_mask, rh_mask) = single(rh_cov_mat);
final_PMM = blkdiag(lh_full_cov_mat, rh_full_cov_mat);
save(fullfile(output_dir, 'premultiplied_matrix_single.mat'), 'final_PMM', 'dim', '-v7.3');
end