function CBIG_MMLDA_vector2nifti(vector, mask, out_name)
% CBIG_MMLDA_vector2nifti(vector, mask, out_name)
%
% Convert vector to a nifti volume by filling the voxels in the mask of value 1 with values in the vector.
%
% Input:
%   - vector    : 1 x M vector. M is number of voxels.
%   - mask      : nifti volume. Number of 1s in mask should be M.
%   - out_name  : output name of nifti volume with '.nii.gz' suffix
%
% Example:
%   CBIG_MMLDA_vector2nifti(vector, 'mask', '~/factor1.nii.gz')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

mri = MRIread(mask);
gm_mask = mri.vol;
mask_size = size(gm_mask);
gm_mask(gm_mask~=0) = 1;
gm_mask_1d = reshape(gm_mask, [1 mask_size(1)*mask_size(2)*mask_size(3)]);

vol_row_1d = gm_mask_1d;
vol_row_1d(vol_row_1d==1) = vector;
vol_row_3d = reshape(vol_row_1d, [mask_size(1) mask_size(2) mask_size(3)]);

mri.vol = vol_row_3d;
MRIwrite(mri, out_name);