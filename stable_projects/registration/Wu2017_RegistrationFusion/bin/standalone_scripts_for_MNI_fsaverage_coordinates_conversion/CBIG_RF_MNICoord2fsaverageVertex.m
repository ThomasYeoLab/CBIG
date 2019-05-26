function [vertices, fs_coords] = CBIG_RF_MNICoord2fsaverageVertex(mni_coords, surf_mesh)
% [vertices, fs_coords] = CBIG_RF_MNICoord2fsaverageVertex(mni_coords, surf_mesh)
%
% This function takes in RAS coordinates in MNI152 space, to convert to the corresponding 
% vertex number and coordinates in fsavearge surface space. 
%
% Note that this conversion uses the fsaverage-to-MNI152 mapping. Therefore the results 
% may not be consistent with using the final projection scripts (or the standalone scripts)
% to project MNI152 data to fsaverage surface.
%
% Input:
%     - mni_coords:
%                  3xN matrix containing the N set of RAS coordinates in MNI152 space to 
%                  convert
%     - surf_mesh :
%                  type of surface mesh to use for fsaverage vertex coordinates ('white', 
%                  'pial', 'sphere', or 'inflated')
%                  (default: 'white')
%
% Output:
%     - vertices :
%                 1xN array of the N corresponding vertex numbers on fsaverage surface
%     - fs_coords:
%                 3xN matrix containing the coordinates of the N vertices in the selected surface mesh
%
% Example:
% 1) vertices = CBIG_RF_MNICoord2fsaverageVertex([60; 0; 10])
% This command finds the corresponding vertex numbers on fsaverage surface for the voxel at (60, 0, 10) 
% in MNI152 space. This vertices output should be 9092.
%
% 2) vertices = CBIG_RF_MNICoord2fsaverageVertex([0; 0; 0])
% This command finds the corresponding vertex numbers on fsaverage surface for the voxel at (0, 0, 0) 
% in MNI152 space. As this voxel is not in the cortex, a warning would be printed and the output would be nan.
%
% Note that you should be in the same directory as this script to run the examples.
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Function usage
if nargin < 1
    disp('usage: [vertices, fs_coords] = CBIG_RF_MNICoord2fsaverageVertex(mni_coords, surf_mesh)');
    return
end

% Default parameter
if nargin < 2
    surf_mesh = 'white';
end

% Load mappings
dir_uti = fileparts(fileparts(mfilename('fullpath')));
map = fullfile(dir_uti, 'final_warps_FS5.3', 'allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.vertex.mat');
load(map, 'lh_vertex', 'rh_vertex');

% Mask out non-cortical areas
mask = MRIread(fullfile(dir_uti, 'liberal_cortex_masks_FS5.3', 'FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz'));
lh_vertex(mask.vol==0) = 0;
rh_vertex(mask.vol==0) = 0;

% Convert input RAS coordinates to matrix coordiantes
mni_vox = convertRas2Vox(mni_coords, mask.vox2ras);
mni_mat = [mni_vox(2, :)+1; mni_vox(1, :)+1; mni_vox(3, :)+1];

% Read surface mesh
lh_mesh = read_surf(fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage', 'surf', ['lh.' surf_mesh]));
rh_mesh = read_surf(fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage', 'surf', ['rh.' surf_mesh]));

% Get the corresponding vertices and coordinates
n = size(mni_coords, 2);
vertices = zeros(1, n);
fs_coords = zeros(3, n);
for i = 1:n
    lh_corr = lh_vertex(mni_mat(1, i), mni_mat(2, i), mni_mat(3, i));
    rh_corr = rh_vertex(mni_mat(1, i), mni_mat(2, i), mni_mat(3, i));
    if lh_corr ~= 0 % grab vertex from left hemisphere
        vertices(i) = lh_corr;
        fs_coords(:, i) = lh_mesh(vertices(i), :);
    elseif rh_corr ~= 0 % grab vertex from right hemisphere
        vertices(i) = rh_corr;
        fs_coords(:, i) = rh_mesh(vertices(i), :);
    else % vertex is outside loose cortical mask
        fprintf('Coordinate set %d does not fall inside our liberal cortical mask. ', i);
        fprintf('No reasonable correspondence could be established. ');
        fprintf('The output for set %d will be marked as NaN. \n', i);
        vertices(i) = nan;
        fs_coords(:, i) = [nan nan nan]';
    end       
end

end

function vox = convertRas2Vox(ras, vox2ras)

vox = inv(vox2ras(1:3, 1:3)) * (ras - repmat(vox2ras(1:3, 4), 1, size(ras, 2)));

end
