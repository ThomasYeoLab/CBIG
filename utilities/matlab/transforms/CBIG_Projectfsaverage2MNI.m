function output = CBIG_Projectfsaverage2MNI(lh_surf, rh_surf, MNI_mask, index_volume, delimiter)

% Usage: output = CBIG_Projectfsaverage2MNI(lh_surf, rh_surf, MNI_mask, index_volume, delimiter)
%
% Warning! This function is obsolete. Please use CBIG_Projectfsaverage2MNI_Ants instead.
%
% Function projects surface data to the volume
% Also see reverse transform CBIG_ProjectMNI2fsaverage2
% 
% 
% ------------------------------------------------------------
% EXAMPLE USAGE: Project cortical gyral labels to MNI template
% ------------------------------------------------------------
% >> lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'aparc.a2009s.annot');
% >> output = CBIG_Projectfsaverage2MNI(lh_avg_mesh.MARS_label', rh_avg_mesh.MARS_label');
% >> MRIwrite(output, 'test.nii.gz');
% 
%
% -----------
% DEFINITIONS
% -----------
% lh_surf      = data matrix from left hemi : num_vertices x data_dimension
% rh_surf      = data matrix from right hemi: num_vertices x data_dimension   
% output_file  = volumetric output filename
% MNI_mask     = filename of binary mask in MNI space 
%                (default = $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Buckner2011/coord_vol2surf/MNI_cortex_estimate.150.nii.gz)
%                If set to 'NONE' means do not specify a mask 
% index_volume = filename of volume that specifies for each voxel in volumetric space, corresponding vertex number in surface space. 
%                Positive indices correspond to left hemi vertices. Negative indices correspond to right hemi vertices. 
%                e.g., value of 23 corresponds to vertex 23 on left hemi. value of -23 corresponds to vertex 23 of right hemi.
%                (default = $CBIG_CODE_DIR/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Buckner2011/coord_vol2surf/1000sub.FSL_MNI152.1mm.full_vertex_map.500.nii.gz')
% delimiter    = values outside MNI mask is set to delimiter value (default = 0)
% 
% The default MNI_mask and index_volume corresponds to MNI template 
% $CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz, which is 
% the FSL MNI 1mm template "freesurfer conformed" to 256 x 256 x 256 resolution.
%
%
% ------------------------
% CREATION OF INDEX_VOLUME
% ------------------------
% The default index_volume is computed by running 1000 subjects and FSL MNI152 template through recon-all. 
% Warps between fsaverage <-> each subject native anatomical space <->
% freesurfer nonlinear volumetric space <-> MNI152 template are then
% averaged across 1000 subjects. Since only cortical voxels have valid surface correspondence, 
% we also keep track of the number of subjects, whose cortical surface maps to a single voxel.
% For voxels which are mapped onto by less than 500 subjects, we assign them a cortical membership based on the closest voxel which does. 
%
%
% ------------------------
% CREATION OF MNI MASK
% ------------------------
% The default MNI_mask is a loose cortical mask. Running the MNI template through the recon-all pipeline 
% gives a cortex.nii.gz mask which severe underestimates cortical voxels. We produce MNI_cortex_estimate.150.nii.gz as follows
% 
%   1) Create intermediate MNI cortical mask, where a voxel is considered a cortical voxel if the cortex 
%      of at least 150 subjects maps to the voxel OR if recon-all of MNI template decides the voxel is a
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
%   7) Visual inspection to see if it looks ok --> MNI_cortex_estimate.150.nii.gz
%
% ----------
% References
% ----------
%     1) Yeo BTT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari L, Hollinshead M, Roffman JL, Smoller JW, Z�llei L, Polimeni JM, Fischl B, Liu H, Buckner RL. 
%        The organization of the human cerebral cortex revealed by intrinsic functional connectivity. 
%        J Neurophysiology, 106(3):1125?1165, 2011
%
%     2) Buckner RL, Krienen FM, Castellanos A, Diaz JC, Yeo BTT. 
%        The organization of the human cerebellum revealed by intrinsic functional connectivity. 
%        J Neurophysiology, 106(5):2322-2345, 2011
% 
%     3) Choi EY, Yeo BTT, Buckner RL.
%        The organization of the human striatum revealed by intrinsic functional connectivity. 
%        J Neurophysiology, 108(8):2242-2263, 2012 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

warning('This function is obsolete. Please use CBIG_Projectfsaverage2MNI_Ants instead.');

if(nargin < 3)
    % loose mask
    MNI_mask = fullfile(getenv('CBIG_CODE_DIR'), '/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Buckner2011/coord_vol2surf/', 'MNI_cortex_estimate.150.nii.gz');
end

if(nargin < 4)
    index_volume = fullfile(getenv('CBIG_CODE_DIR'), '/utilities/matlab/transforms/CorrespondenceFreeSurferVolSurfSpace_Buckner2011/coord_vol2surf/', '1000sub.FSL_MNI152.1mm.full_vertex_map.500.nii.gz');
end

if(nargin < 5)
    delimiter = 0;
else
    if(ischar(delimiter))
        delimiter = str2num(delimiter);
    end
end

% read index volume
index = MRIread(index_volume);
pos_index = find(index.vol > 0);
neg_index = find(index.vol < 0);

% Check index volume
if(sum(index.vol == 0) > 0)
   error('There are index with 0 values'); 
end

if((max(abs(index.vol(:))) > size(lh_surf, 1)) || (max(abs(index.vol(:))) > size(rh_surf, 1)))
    error('Index volume is indexing into surface that is higher resolution than surface data');
end

% read mask
if(~strcmp(MNI_mask, 'NONE'))
    mask = MRIread(MNI_mask);
    
    if((size(mask.vol, 1) ~= size(index.vol, 1)) || (size(mask.vol, 2) ~= size(index.vol, 2)) || (size(mask.vol, 3) ~= size(index.vol, 3))) 
       error('mask not the same size as index volume'); 
    end
    
    if(max(abs(mask.vox2ras(:) - index.vox2ras(:))) > 1e-5)
       warning('mask and index volume has different header information. Are you sure they are in the same space?');
    end
end

% Perform projection
output = index;
output.vol = zeros([size(index.vol) size(lh_surf, 2)]);
dummy = zeros(size(index.vol));
for i = 1:size(lh_surf, 2)
    dummy(pos_index) = lh_surf(index.vol(pos_index), i);
    dummy(neg_index) = rh_surf(abs(index.vol(neg_index)), i);
    
    if(~strcmp(MNI_mask, 'NONE'))
        dummy(mask.vol == 0) = delimiter;
    end
    output.vol(:, :, :, i) = dummy;
end







