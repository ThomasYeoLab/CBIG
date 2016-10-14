function [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage(data, average)

% This function is specific for Buckner2011 paper, please use more general
% function CBIG_ProjectMNI2fsaverage2.m
%
% [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage(data, average)
% assumes data is read from MRIread
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

 

data_dim = size(data.vol, 4);
for hemis = {'lh', 'rh'}
    hemi = hemis{1};
    load(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Buckner2011', 'coord_surf2vol', [hemi '.1000sub.FSL_MNI.ras.mat']));

    avg_mesh = CBIG_ReadNCAvgMesh(hemi, average, 'inflated', 'cortex');
    num_vertices = size(avg_mesh.vertices, 2);
    proj_data = zeros(data_dim, num_vertices);
    
    % Convert RAS correspondence to voxel coordinates and matrix coordinates
    vox_coor = CBIG_ConvertRas2Vox(ras(:, 1:num_vertices), data.vox2ras);
    mat_coor = [vox_coor(2, :)+1 ; vox_coor(1, :)+1; vox_coor(3, :)+1];
    
    % Convert matrix coordinates correspondence to index
    mat_coor = round(mat_coor);
    
    for i = 1:data_dim
        mat_index = sub2ind(size(data.vol), mat_coor(1, :), mat_coor(2, :), mat_coor(3, :), i*ones(1, num_vertices));
        proj_data(i, :) = data.vol(mat_index);
    end
    
    disp(['# vertices: ' num2str(num_vertices)]);
    disp(['# unique voxels: ' num2str(length(unique(mat_index)))]);
    
    if(strcmp(hemi, 'lh'))
        lh_proj_data = proj_data;
    else
        rh_proj_data = proj_data;
    end
end
