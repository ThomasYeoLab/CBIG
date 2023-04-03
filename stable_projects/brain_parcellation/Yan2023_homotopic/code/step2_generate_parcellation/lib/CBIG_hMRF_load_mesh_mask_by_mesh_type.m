function [lh_avg_mesh, rh_avg_mesh, lh_mask, rh_mask] = CBIG_hMRF_load_mesh_mask_by_mesh_type(mesh_type)
% [lh_avg_mesh, rh_avg_mesh, lh_mask, rh_mask] = CBIG_hMRF_load_mesh_mask_by_mesh_type(mesh_type)
%
% This function loads the surface mesh and the binary vector indicating cortex vs medial wall given a surface mesh.

% For the notations below:
% N = no of vertices per hemisphere; 

% Input
%   - mesh_type: (string)
%     The string indicating the desired mesh type.
% Output
%   - lh_avg_mesh: (struct)
%     The structure containing mesh information for the left hemisphere.
%   - rh_avg_mesh: (struct)
%     The structure containing mesh information for the right hemisphere.
%   - lh_mask: (matrix)
%     A Nx1 binary array of length N (total vertices of the given surface mesh) indicating whether a vertex belongs 
%     to cortex or the medial wall on the left hemisphere.
%   - rh_mask: (matrix)
%     A Nx1 binary array of length N (total vertices of the given surface mesh) indicating whether a vertex belongs 
%     to cortex or the medial wall on the right hemisphere.
%
% Example
%   - [lh_avg_mesh, rh_avg_mesh, lh_mask, rh_mask] = CBIG_hMRF_load_mesh_mask_by_mesh_type('fsaverage6')
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(contains(mesh_type, 'fsaverage'))
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_type, 'inflated', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_type, 'inflated', 'cortex');
    lh_mask = (lh_avg_mesh.MARS_label == 2)';
    rh_mask = (rh_avg_mesh.MARS_label == 2)';
elseif(strcmp(mesh_type, 'fs_LR_32k'))
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh_type, 'very_inflated','medialwall.annot');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh_type, 'very_inflated','medialwall.annot');
    lh_mask = ~(lh_avg_mesh.MARS_label==1);
    rh_mask = ~(rh_avg_mesh.MARS_label==1);
else
    error('Input mesh not recognized!');
end
end