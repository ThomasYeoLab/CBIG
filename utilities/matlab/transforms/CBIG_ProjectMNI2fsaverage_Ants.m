function [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage_Ants(data, average, interp_type, lh_index_mat, rh_index_mat)

% Usage: [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage_Ants(data, average, interp_type, lh_index_mat, rh_index_mat)
%
% Function projects MNI data to fsaverage
% Also see reverse transform CBIG_Projectfsaverage2MNI_Ants
% 
% ---------------------------------------------------------------------
% EXAMPLE USAGE: Project cortical gyral labels to MNI template and back
% ---------------------------------------------------------------------
% >> lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> CBIG_DrawSurfaceMaps(lh_avg_mesh.MARS_label, rh_avg_mesh.MARS_label, 'fsaverage', 'inflated', 0, 36);
% >> data = CBIG_Projectfsaverage2MNI_Ants(lh_avg_mesh.MARS_label', rh_avg_mesh.MARS_label');
% >> [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage_Ants(data, 'fsaverage', 'nearest'); 
% >> CBIG_DrawSurfaceMaps(lh_proj_data(1, :), rh_proj_data(1, :), 'fsaverage', 'inflated', 0, 36); 
% 
% -----------
% INTPUTS
% -----------
%    - data        : 
%                    assumed to be read using MRIread in MNI152 space (can be 4 dimensional)
%    - average     : 
%                    'fsaverage', 'fsaverage6' or 'fsaverage5'
%    - interp_type : 
%                    'linear' (default), 'nearest', 'spline', 'cubic'
%    - lh_index_mat: 
%                    .mat file containing matrix "ras": a 3 x num_vertices matrix specifying for each left hemi vertex 
%                    the corresponding MNI152 RAS coordinates. 
%                    (default is $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Wu2017/final_warps_FS5.3/lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat)
%    - rh_index_mat: 
%                    .mat file containing matrix "ras": a 3 x num_vertices matrix specifying for each right hemi vertex 
%                    the corresponding MNI152 RAS coordinates. 
%                    (default is $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Wu2017/final_warps_FS5.3/rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat)       
% 
% The default index matrices are 3 x 163842 corresponding to fsaverage
% resolution. Even if you are interested in lower resolution, I suggest you
% project to fsaverage resolution, smooth and downsample.
%
% Due to the way the code works, the data volume just needs to be in MNI152
% space. The data can be 2mm, 3mm and with different headers. That's ok as long as
% their corresponding RAS coordinates are correct. To check, I recommend
% using freeview to load your data volume and the MNI152 template from
% here: $FSL_DIR/data/standard/MNI152_T1_1mm_brain.nii.gz
% 
% -----------
% OUTPUTS
% -----------
%    - lh_surf: 
%               projected results in fsaverage from left hemisphere 
%               size = num_vertices x data_dimension
%    - rh_surf: 
%               projected results in fsaverage from right hemisphere
%               size = num_vertices x data_dimension
%
% --------------------------------------
% CREATION OF LH_INDEX_MAT, RH_INDEX_MAT
% --------------------------------------
% The default index_mats are computed by running 1490 GSP subjects through recon-all.
% Subjects' native anatomical spaces are registered to FSL MNI152 template by ANTs
% Warps between fsaverage <-> each subject native anatomical space <-> MNI152 template are then averaged across the 1490 subjects. 
%     
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%Set up default parameters
if(nargin < 3)
   interp_type = 'linear'; 
end

if(nargin < 4)
    lh_index_mat = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Wu2017', 'final_warps_FS5.3', 'lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
end

if(nargin < 5)
    rh_index_mat = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Wu2017', 'final_warps_FS5.3', 'rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat');
end

%Set up grid index in data.vol space
[x1, x2, x3] = ndgrid(1:size(data.vol, 1), 1:size(data.vol, 2), 1:size(data.vol, 3));

%Loop through each hemisphere
data_dim = size(data.vol, 4);
for hemis = {'lh', 'rh'}
    hemi = hemis{1};
    
    %Load corresponding hemisphere's warp
    if(strcmp(hemi, 'lh'))
        load(lh_index_mat);
    else
        load(rh_index_mat);
    end

    %Read number of vertices of fsaverage surface
    avg_mesh = CBIG_ReadNCAvgMesh(hemi, average, 'inflated', 'cortex');
    num_vertices = size(avg_mesh.vertices, 2);
    proj_data = zeros(data_dim, num_vertices);
    
    % Convert RAS correspondence to voxel coordinates and matrix coordinates
    vox_coor = CBIG_ConvertRas2Vox(ras(:, 1:num_vertices), data.vox2ras);
    mat_coor = [vox_coor(2, :)+1 ; vox_coor(1, :)+1; vox_coor(3, :)+1];
    
    %Project volume to fsaverage
    %Loop through each frame for 4D volume
    for i = 1:data_dim
        proj_data(i, :) = interpn(x1, x2, x3, squeeze(data.vol(:, :, :, i)), mat_coor(1, :)', mat_coor(2, :)', mat_coor(3, :)', interp_type)';
    end
    disp(['# vertices: ' num2str(num_vertices)]);
    
    %Assign projected data to output
    if(strcmp(hemi, 'lh'))
        lh_proj_data = proj_data;
    else
        rh_proj_data = proj_data;
    end
end
