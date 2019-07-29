function [lh_label_fslr, rh_label_fslr, lh_label_fslr_164k, rh_label_fslr_164k] = ...
  CBIG_project_from_fsaverage_to_fslr(lh_label, rh_label, abs_path_to_output_folder, type_of_smoothing)

% [lh_label_fslr, rh_label_fslr, lh_label_fslr_164k, rh_label_fslr_164k] =
%  CBIG_project_from_fsaverage_to_fslr(lh_label, rh_label, abs_path_to_output_folder, type_of_smoothing)
%
% This function projects surface labels (parcellation labels, activation
% patterns...) from fsaverage5/fsaverage6/fsaverage to fslr_32k/fslr_164k
% Input:
%      - lh_label, rh_label:
%        vector of surface labels on the left and right hemishere respectively in
%        fsaverage5/fsaverage6/fsaverage.
%        Each variable is a Nx1 vectore, where:
%        N = 40962 for fsaverage6
%        N = 10242 for fsaverage5
%        N = 163842 for fsaverage
%
%      - abs_path_to_output_folder:
%        absolute path to an output folder.
%
%      - type_of_smoothing:
%        'METRIC_NEAREST_NODE' for integer label data.
%        'METRIC_AVERAGE_TILE' for other types of label data.
%        If not given, the default is 'METRIC_NEAREST_NODE'.
% Output:
%      - lh_label_fslr,rh_label_fslr:
%        projected surface labels of the left and right hemisphere respectively
%        on fs_LR_32k surface.
%        Each variable is a Nx1 vector, where N = 32492.
%
%      - lh_label_fslr_164k,rh_label_fslr_164k:
%        projected surface labels of the left and right hemisphere respectively
%        on fs_LR_164k surface.
%        Each variable is a Nx1 vector, where N = 163842.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(lh_label,2) ~= 1
    error('Input argument ''lh_label'' should be a column vector');
end
if size(rh_label,2) ~= 1
    error('Input argument ''rh_label'' should be a column vector');
end

if (nargin<3) % if no output folder is provided
    error('Please provide the absolute path to an output folder')
end
if (nargin<4)
    type_of_smoothing='METRIC_NEAREST_NODE'; %default is for integer label data
    warning('No type of smoothing is provided, using default as METRIC_NEAREST_NODE')
end

% unclear why, but this might cause problems otherwise
if (size(lh_label,1) ~= 1);
    lh_label=lh_label';
end

if (size(rh_label,1) ~= 1);
    rh_label=rh_label';
end

if (max(size(lh_label)) == 40962) % if input is fsaverage6 data, we upsample from fsaverag6 to fsaverage
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');
    rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'white', 'cortex');
    lh_mesh6 = CBIG_ReadNCAvgMesh('lh', 'fsaverage6', 'white', 'cortex');
    rh_mesh6 = CBIG_ReadNCAvgMesh('rh', 'fsaverage6', 'white', 'cortex');
    
    lh_labels7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices,lh_mesh6,lh_label);
    rh_labels7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices,rh_mesh6,rh_label);
elseif (max(size(lh_label)) == 10242) % if input is fsaverage5 data, we upsample from fsaverag5 to fsaverage
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');
    rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'white', 'cortex');
    lh_mesh5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'white', 'cortex');
    rh_mesh5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'white', 'cortex');
    
    lh_labels7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices,lh_mesh5,lh_label);
    rh_labels7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices,rh_mesh5,rh_label);
elseif (max(size(lh_label)) == 163842) % if input is fsaverage data, no upsampling is needed
    lh_labels7=lh_label;
    rh_labels7=rh_label;
else
    error('Can only perform projection on fsaverage5/fsaverage6/fsaverage input data')
end

% start processing

if (strcmp(type_of_smoothing, 'METRIC_NEAREST_NODE'))
    lh_mri.vol = int16(lh_labels7);
    rh_mri.vol = int16(rh_labels7);
    MRIwrite(lh_mri, fullfile(abs_path_to_output_folder, 'label_lh_borders.mgh'), 'int');
    MRIwrite(rh_mri, fullfile(abs_path_to_output_folder, 'label_rh_borders.mgh'), 'int');
elseif (strcmp(type_of_smoothing, 'METRIC_AVERAGE_TILE')) % assuming double data
    warning('Assuming input data is of double type')
    lh_mri.vol = double(lh_labels7);
    rh_mri.vol = double(rh_labels7);
    MRIwrite(lh_mri, fullfile(abs_path_to_output_folder, 'label_lh_borders.mgh'), 'double');
    MRIwrite(rh_mri, fullfile(abs_path_to_output_folder, 'label_rh_borders.mgh'), 'double');
else
    error('Unknown type of smoothing')
end

% status and cmdout are saved for debugging
SCRIPT_PATH = fullfile(getenv('CBIG_CODE_DIR') , ...
  'utilities/matlab/fslr_matlab/bash_scripts/CBIG_transfer_border_to_fs_lr_32k.sh');
[status,cmdout] = system([SCRIPT_PATH ' ' abs_path_to_output_folder ' label ' type_of_smoothing]);

rh_gifti = gifti(fullfile(abs_path_to_output_folder, 'label.R.32k_fs_LR.func.gii'));
lh_gifti = gifti(fullfile(abs_path_to_output_folder, 'label.L.32k_fs_LR.func.gii'));
lh_label_fslr = lh_gifti.cdata;
rh_label_fslr = rh_gifti.cdata;

rh_gifti = gifti(fullfile(abs_path_to_output_folder, 'label.R.fs_LR.func.gii'));
lh_gifti = gifti(fullfile(abs_path_to_output_folder, 'label.L.fs_LR.func.gii'));
lh_label_fslr_164k = lh_gifti.cdata;
rh_label_fslr_164k = rh_gifti.cdata;

end
