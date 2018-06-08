function output = CBIG_Projectfsaverage2MNI_Ants(lh_surf, rh_surf, interp_type, MNI_mask, index_map, delimiter)

% Usage: output = CBIG_Projectfsaverage2MNI_Ants(lh_surf, rh_surf, interp_type, MNI_mask, index_map, delimiter)
%
% Function projects surface data to the volume
% Also see reverse transform CBIG_ProjectMNI2fsaverage_Ants
% 
% ------------------------------------------------------------
% EXAMPLE USAGE: Project cortical gyral labels to MNI template
% ------------------------------------------------------------
% >> lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> output = CBIG_Projectfsaverage2MNI_Ants(lh_avg_mesh.MARS_label', rh_avg_mesh.MARS_label');
% >> MRIwrite(output, 'test.nii.gz');
%
%
% -----------
% INPUTS
% -----------
%    - lh_surf    : 
%                   data matrix from left hemi : num_vertices x data_dimension
%    - rh_surf    : 
%                   data matrix from right hemi: num_vertices x data_dimension
%    - interp_type: 
%                   'linear' (default) or 'nearest'
%    - MNI_mask   : 
%                   filename of binary mask in MNI space 
%                   (default = $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Wu2017/liberal_corte_masks_FS5.3/FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz)
%                   If set to 'NONE' means do not specify a mask 
%    - index_map  : 
%                   filename of mapping that specifies for each voxel in volumetric space, corresponding spherical coordinates in surface space. 
%                   (default = $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Wu2017/final_warps_FS5.3/allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.prop.mat')
%    - delimiter  : 
%                   values outside MNI mask is set to delimiter value (default = 0)
%                   note that delimiter should be passed in as char, e.g. '1000'
% 
% The default MNI_mask and index_volume corresponds to MNI template 
% $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz, which is 
% the FSL MNI 1mm template "freesurfer conformed" to 256 x 256 x 256 resolution.
%
% -----------
% OUTPUTS
% -----------
%    - output    : 
%                  projected data in MNI152 space 
%                  size = 256 x 256 x 256 x data_dimension
%
% ------------------------
% CREATION OF INDEX_VOLUME
% ------------------------
% The default index_volume is computed by running 1490 GSP subjects through recon-all. 
% Subjects' native anatomical spaces are registered to FSL MNI152 template by ANTs
% Warps between fsaverage <-> each subject native anatomical space <-> MNI152 template 
% are then averaged across 1490 subjects. Since only cortical voxels have valid surface correspondence, 
% we also keep track of the number of subjects, whose cortical surface maps to a single voxel.
% For voxels which are mapped onto by less than 15% of the subjects, we assign them a cortical membership based on the closest voxel which does. 
%
%
% ------------------------
% CREATION OF MNI MASK
% ------------------------
% The default MNI_mask is a loose cortical mask. Running the MNI template through the recon-all pipeline 
% gives a cortex.nii.gz mask which severe underestimates cortical voxels. We produce FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz as follows
% 
%   1) Create intermediate MNI cortical mask, where a voxel is considered a cortical voxel if the cortex 
%      of at least 15% of all 1490 subjects maps to the voxel OR if recon-all of MNI template decides the voxel is a
%      cortical voxel. 
% 
%   2) Smooth intermediate mask: smooth3(double(mask), 'box', 5)
%
%   3) Threshold at 0.5
% 
%   4) Use aparc+aseg of MNI template to mask out voxels that are for sure non cortex, i.e, 
%        aparc_aseg = MRIread('$CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS/mri/aparc+aseg.mgz');
%        aparc_aseg = ((aparc_aseg.vol < 1000) & (aparc_aseg.vol > 0) & (aparc_aseg.vol ~= 41) & (aparc_aseg.vol ~= 2));
%        mask(aparc_aseg) = 0;
% 
%   5) Fill holes: imfill(mask, 'holes')
%
%   6) Remove islands using bwlabeln and removing islands less than 100
%
%   7) Visual inspection to see if it looks ok --> FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz
% 
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%Set up default parameters
if(nargin < 3)
   interp_type = 'linear'; 
end

if(nargin < 4)
    % loose mask
    MNI_mask = fullfile(getenv('CBIG_CODE_DIR'), '/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Wu2017/liberal_cortex_masks_FS5.3/', 'FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz');
end

if(nargin < 5)
    index_map = fullfile(getenv('CBIG_CODE_DIR'), '/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Wu2017/final_warps_FS5.3/', 'allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.prop.mat');
end

if(nargin < 6)
    delimiter = 0; %set voxels outside loose mask to 0 by default
else
    if(ischar(delimiter)) %otherwise user should set a char value
        delimiter = str2num(delimiter);
    else
        error('Invalid delimiter.');
    end
end

% read index mapping
load(index_map);
%Voxels with (0, 0, 0) mapping will be masked out
lh_mask = double(sum(abs(lh_coord))~=0);
rh_mask = double(sum(abs(rh_coord))~=0);

%Load fsaverage spherical cortex mesh
if size(lh_surf, 1) == 163842
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex');
else
    error('Invalid number of vertices. If your data is in fsaverage5/fsaverage6, please upsample to fsaverage first.');
end
    
% read mask
if(~strcmp(MNI_mask, 'NONE'))
    mask = MRIread(MNI_mask);
end

%Project source surface to volume space by the chosen interpolation
output = mask;
output.vol = zeros([size(mask.vol) size(lh_surf, 2)]);
for i = 1:size(lh_surf, 2)
    switch interp_type
        case 'nearest'
            %Project lh
            lh_vertex = zeros(1, size(lh_coord, 2));
            lh_vertex(lh_mask~=0) = MARS_findNV_kdTree(single(lh_coord(:, lh_mask~=0)), lh_avg_mesh.vertices);
            lh_projected = zeros(1, size(lh_coord, 2));
            lh_projected(lh_mask~=0) = interpn(1:size(lh_surf, 1), lh_surf(:, i)', lh_vertex(lh_mask~=0), 'nearest');
            
            %Project rh
            rh_vertex = zeros(1, size(rh_coord, 2));
            rh_vertex(rh_mask~=0) = MARS_findNV_kdTree(single(rh_coord(:, rh_mask~=0)), rh_avg_mesh.vertices);
            rh_projected = zeros(1, size(lh_coord, 2));
            rh_projected(rh_mask~=0) = interpn(1:size(rh_surf, 1), rh_surf(:, i)', rh_vertex(rh_mask~=0), 'nearest');
            
        case 'linear'
            %Project lh
            lh_projected = zeros(1, size(lh_coord, 2));
            lh_projected(lh_mask~=0) = MARS_linearInterpolate_kdTree(single(lh_coord(:, lh_mask~=0)), lh_avg_mesh, single(lh_surf(:, i)'));
            
            %Project rh
            rh_projected = zeros(1, size(rh_coord, 2));
            rh_projected(rh_mask~=0) = MARS_linearInterpolate_kdTree(single(rh_coord(:, rh_mask~=0)), rh_avg_mesh, single(rh_surf(:, i)'));
        otherwise
            disp('Invalid interpolation option. Use either nearest or linear.');
    end
    
    %Mask out non-cortical areas using the mask provided
    lh_error = sum(abs(double(lh_projected~=0) - mask.vol(:)'), 2);
    rh_error = sum(abs(double(rh_projected~=0) - mask.vol(:)'), 2);
    if(~strcmp(MNI_mask, 'NONE'))
        lh_projected(mask.vol(:)==0) = delimiter;
        rh_projected(mask.vol(:)==0) = delimiter;
    end
    disp(['Total error for lh = ' num2str(lh_error) ', rh = ' num2str(rh_error)]);
    
    %Combine results
    output.vol(:, :, :, i) = reshape(lh_projected + rh_projected, size(mask.vol));
end
