function [lh_label, rh_label] = CBIG_project_from_fslr_to_fsaverage(lh_label_fslr, rh_label_fslr, abs_path_to_output_folder)

% [lh_label, rh_label] = CBIG_project_from_fslr_to_fsaverage(lh_label_fslr, rh_label_fslr, abs_path_to_output_folder)
%
% This function projects labels (parcellation labels, activation 
% patterns...) on a surface in fslr_32k to fsaverage.
%
% Input:
%      - lh_label_fslr, rh_label_fslr:
%        surface labels on the left and right hemisphere respectively in
%        fslr_32k.
%        Each variable is a Nx1 vector, where N = 32,492
%
%      - abs_path_to_output_folder:
%         absolute path to an output folder.
%         Two '.label.gii' files (one for each hemisphere) of labels on
%         fslr_32k surface will be saved to this folder.
%         Two 'func.gii' files (one for each hemisphere) of labels on
%         fsaverage surface will be saved to this folder.
% Output:
%      - lh_label,rh_label:
%        projected surface labels of the left and right hemisphere respectively
%        on fsaverage surface.
%        Each variable is a Nx1 vector, where N = 163842
%
% Example:
% [lh_label, rh_label] = CBIG_project_from_fslr_to_fsaverage(lh_label_fslr, rh_label_fslr, '/data/outputs/fsaverage')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(lh_label_fslr,2) ~= 1
    error('Input argument ''lh_label_fslr'' should be a column vector');
end
if size(rh_label_fslr,2) ~= 1
    error('Input argument ''rh_label_fslr'' should be a column vector');
end

if (nargin<3) % if no output folder is provided
    error('Please provide the absolute path to an output folder')
end

lh_mri.vol = int16(lh_label_fslr);
rh_mri.vol = int16(rh_label_fslr);

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

% get gifti struct
CIFTI_FUNC_DIR = fullfile(CBIG_CODE_DIR, 'data/templates/surface/fs_LR_32k/example_func');
target_lh = gifti(fullfile(CIFTI_FUNC_DIR, 'Parcels_L.func.gii')); % get cifti func file structure
target_rh = gifti(fullfile(CIFTI_FUNC_DIR, 'Parcels_L.func.gii'));

target_lh.cdata = lh_label_fslr;
target_rh.cdata = rh_label_fslr;

save(target_lh, fullfile(abs_path_to_output_folder, 'label.L.32k_fs_LR.label.gii','Base64Binary'));
save(target_rh, fullfile(abs_path_to_output_folder, 'label.R.32k_fs_LR.label.gii','Base64Binary'));

TRANSFORM_SCRIPT = fullfile(CBIG_CODE_DIR, '/utilities/matlab/fslr_matlab/bash_scripts/CBIG_transfer_border_to_fsaverage_164k.sh');
system([TRANSFORM_SCRIPT ' ' abs_path_to_output_folder ' label']);

rh_gifti = gifti(fullfile(abs_path_to_output_folder, 'label.R.fs_average.func.gii'));
lh_gifti = gifti(fullfile(abs_path_to_output_folder, 'label.L.fs_average.func.gii'));
lh_label = lh_gifti.cdata;
rh_label = rh_gifti.cdata;
end
