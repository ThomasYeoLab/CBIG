function vox = CBIG_ConvertRas2Vox(ras, vox2ras)

% Convert RAS coordinate to Voxel coordinate.
% 
%     vox = CBIG_ConvertRas2Vox(ras, vox2ras)
%     Input:
%         ras        : RAS coordinate. 3 x N matrix, where N is num of points.
%         vox2ras    : vox2ras matrix. 4 x 4 matrix, included in nifti file
%     Output:
%         vox        : voxel coordinate. 3 x N matrix, where N is num of points.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


vox = inv(vox2ras(1:3, 1:3)) * (ras - repmat(vox2ras(1:3, 4), 1, size(ras, 2)));


