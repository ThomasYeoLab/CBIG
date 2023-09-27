function [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage2(data, average, interp_type, lh_index_mat, rh_index_mat)

% Usage: [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage2(data, average, interp_type, lh_index_mat, rh_index_mat)
%
% Warning! This function is obsolete. Please use CBIG_ProjectMNI2fsaverage_Ants instead.
%
% Function projects MNI data to fsaverage
% Also see reverse transform CBIG_Projectfsaverage2MNI
% 
%
% ---------------------------------------------------------------------
% EXAMPLE USAGE: Project cortical gyral labels to MNI template and back
% ---------------------------------------------------------------------
% >> lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> CBIG_DrawSurfaceMaps(lh_avg_mesh.MARS_label, rh_avg_mesh.MARS_label, 'fsaverage', 'inflated', 0, 36);
% >> data = CBIG_Projectfsaverage2MNI(lh_avg_mesh.MARS_label', rh_avg_mesh.MARS_label');
% >> [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage2(data, 'fsaverage', 'nearest'); 
% >> CBIG_DrawSurfaceMaps(lh_proj_data(1, :), rh_proj_data(1, :), 'fsaverage', 'inflated', 0, 36);
% 
% Notice the transformation back and forth does not give back the same map.
% In particular, notice the small holes in the return projection? This is because a small number of surface vertices grabs 
% MNI voxels (specified by the default lh_index_mat and rh_index_mat) that turns out to be outside the cortex
% mask specified by MNI_cortex_estimate.150.nii.gz.
% Maybe in a future version, I should expand MNI_cortex_estimate.150.nii.gz to include MNI voxels that the
% surface vertices are mapping to. For now, we can avoid this issue if we
% add the 'NONE' option to CBIG_Projectfsaverage2MNI:
%
% >> lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'aparc.annot');
% >> rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'aparc.annot');
% >> CBIG_DrawSurfaceMaps(lh_avg_mesh.MARS_label, rh_avg_mesh.MARS_label, 'fsaverage', 'inflated', 0, 36));
% >> data = CBIG_Projectfsaverage2MNI(lh_avg_mesh.MARS_label', rh_avg_mesh.MARS_label', 'NONE');
% >> [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage2(data, 'fsaverage', 'nearest'); 
% >> CBIG_DrawSurfaceMaps(lh_proj_data(1, :), rh_proj_data(1, :), 'fsaverage', 'inflated', 0, 36); 
%
% Essentially the 'NONE' option tells CBIG_Projectfsaverage2MNI not to set
% voxels outside the MNI_cortex_estimate.150.nii.gz mask to 0. 
% 
% 
% -----------
% DEFINITIONS
% -----------
% data         = assumed to be read using MRIread in MNI152 space (can be 4 dimensional)
% average      = 'fsaverage', 'fsaverage6' or 'fsaverage5'
% interp_type  = 'linear' (default), 'nearest', 'spline', 'cubic'
% lh_index_mat = .mat file containing matrix "ras" a 3 x num_vertices matrix specifying for each left hemi vertex 
%                the corresponding MNI152 RAS coordinates. 
%                (default is $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Buckner2011/coord_surf2vol/lh.1000sub.FSL_MNI.ras.mat)
% rh_index_mat = .mat file containing matrix "ras" a 3 x num_vertices matrix specifying for each right hemi vertex 
%                the corresponding MNI152 RAS coordinates. 
%                (default is $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Buckner2011/coord_surf2vol/rh.1000sub.FSL_MNI.ras.mat)       
% 
% The default index matrices are 3 x 163842 corresponding to fsaverage
% resolution. Even if you are interested in lower resolution, I suggest you
% project to fsaverage resolution, smooth and downsample.
%
% Due to the way the code works, the data volume just needs to be in MNI152
% space. The data can be 2mm, 3mm and with different headers. That's ok as long as
% their corresponding RAS coordinates are correct. To check, I recommend
% using freeview to load your data volume and the MNI152 template from
% here: $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz
% 
%
% --------------------------------------
% CREATION OF LH_INDEX_MAT, RH_INDEX_MAT
% --------------------------------------
% The default index_mats are computed by running 1000 subjects and FSL MNI152 template through recon-all. 
% Warps between fsaverage <-> each subject native anatomical space <->
% freesurfer nonlinear volumetric space <-> MNI152 template are then
% averaged across 1000 subjects. 
%
%
% ----------
% References
% ----------
%     1) Yeo et al., The organization of the human cerebral cortex estimated by intrinsic functional connectivity.
%         Journal of neurophysiology, 106:1125-1165, 2011
%     
%     2) Buckner et al., The organization of the human cerebellum estimated by intrinsic functional connectivity.
%         Journal of neurophysiology, 106:2322-2345, 2011     
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

warning('This function is obsolete. Please use CBIG_ProjectMNI2fsaverage_Ants instead.');

if(nargin < 3)
   interp_type = 'linear'; 
end

if(nargin < 4)
    lh_index_mat = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Buckner2011', 'coord_surf2vol', 'lh.1000sub.FSL_MNI.ras.mat');
end

if(nargin < 5)
    rh_index_mat = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Buckner2011', 'coord_surf2vol', 'rh.1000sub.FSL_MNI.ras.mat');
end


[x1, x2, x3] = ndgrid(1:size(data.vol, 1), 1:size(data.vol, 2), 1:size(data.vol, 3));

data_dim = size(data.vol, 4);
for hemis = {'lh', 'rh'}
    hemi = hemis{1};
    
    if(strcmp(hemi, 'lh'))
        load(lh_index_mat);
    else
        load(rh_index_mat);
    end

    avg_mesh = CBIG_ReadNCAvgMesh(hemi, average, 'inflated', 'cortex');
    num_vertices = size(avg_mesh.vertices, 2);
    proj_data = zeros(data_dim, num_vertices);
    
    % Convert RAS correspondence to voxel coordinates and matrix coordinates
    vox_coor = CBIG_ConvertRas2Vox(ras(:, 1:num_vertices), data.vox2ras);
    mat_coor = [vox_coor(2, :)+1 ; vox_coor(1, :)+1; vox_coor(3, :)+1];
    
    for i = 1:data_dim
        proj_data(i, :) = interpn(x1, x2, x3, squeeze(data.vol(:, :, :, i)), mat_coor(1, :)', mat_coor(2, :)', mat_coor(3, :)', interp_type)';
    end
    
    disp(['# vertices: ' num2str(num_vertices)]);
    
    if(strcmp(hemi, 'lh'))
        lh_proj_data = proj_data;
    else
        rh_proj_data = proj_data;
    end
end
