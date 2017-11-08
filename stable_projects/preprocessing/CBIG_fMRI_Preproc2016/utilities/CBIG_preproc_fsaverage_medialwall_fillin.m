function CBIG_preproc_fsaverage_medialwall_fillin(hemi, surf_mesh, origin_fs_nii, smooth_fs_nii, output_fs_nii)

%CBIG_preproc_fsaverage_medialwall_fillin(hemi, surf_mesh, origin_fs_nii, smooth_fs_nii)
% This function is used to extract the medial wall data from the original
% surface nifti file, and fill it into the medial wall area of the nifti
% file after smoothing. The reason is FreeSurfer 5.3.0 will set the medial
% wall area to 0 when using mri_surf2surf to perform  smoothing. If we
% downsample the smoothed nifti file to a lower resolution, the medial wall
% after downsampling may not be the area defined by lower resolution
% surface template.
%
% INPUT:                                                                   
%      -hemi:           
%       'lh' or 'rh'                  
%
%      -surf_mesh:      
%       surface mesh of the original surface nifti file, e.g.'fsaverage6'
%
%      -origin_fs_nii:  
%       original surface nifti file before smoothing, e.g. *_fs6.nii.gz
%
%      -smooth_fs_nii: 
%       nifti file after smoothing the original surface nifti file,
%       e.g. *_fs6_sm6.nii.gz.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% Read in the original surf_mesh nifti file and the nifti file after smoothing
origin_surf = MRIread(origin_fs_nii);
smooth_surf = MRIread(smooth_fs_nii);

%reshape the data into a N x T matrix, where N is the number of surface vertices, T is the number of timepoints
origin_data = reshape(origin_surf.vol, size(origin_surf.vol,1)*size(origin_surf.vol,2)*size(origin_surf.vol,3),size(origin_surf.vol,4));
smooth_data = reshape(smooth_surf.vol, size(smooth_surf.vol,1)*size(smooth_surf.vol,2)*size(smooth_surf.vol,3),size(smooth_surf.vol,4));
 
%% Read in surf_mesh cortex label from template and create mask for current hemisphere 
hemi_cortex_info = read_label(surf_mesh,[hemi '.cortex']); 

%cortex label is the first column, label index starts from 0, matlab index need + 1
hemi_cortex_label = hemi_cortex_info(:,1)+1;

%hemi_cortex_mask is a Nx1 vector, 0 indicates cortex, 1 indicates medial wall
hemi_cortex_mask = ones(size(smooth_data,1),1);
hemi_cortex_mask(hemi_cortex_label) = 0;

%% Apply the mask on original surf_mesh nifti file to extract the medial wall data and fill it in the medial wall of smoothed nifti file
smooth_data(hemi_cortex_mask==1,:) = smooth_data(hemi_cortex_mask==1,:) + origin_data(hemi_cortex_mask==1,:);

%%reshape the smoothed data after filling medial wall
smooth_surf.vol = reshape(smooth_data, size(origin_surf.vol,1),size(origin_surf.vol,2),size(origin_surf.vol,3),size(origin_surf.vol,4));

%% Write the smooth data into the smooth_fs_nii file so that the medial wall area will be filled with the data from the origin_fd_nii file
MRIwrite(smooth_surf, output_fs_nii);
