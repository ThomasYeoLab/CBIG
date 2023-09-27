function ras = CBIG_ConvertVox2Ras(vox, vox2ras)

% Convert Voxel coordinate to RAS coordinate.
% 
%     ras = CBIG_ConvertVox2Ras(vox, vox2ras)
%     Input:
%         vox        : voxel coordinate. 3 x N matrix, where N is num of points.
%         vox2ras    : vox2ras matrix. 4 x 4 matrix, included in nifti file
%     Output:
%         ras        : RAS coordinate. 3 x N matrix, where N is num of points.
%
% PS: vox is arranged in Freeview sequence, therefore mat(i, j, k) is
% equivalent to vox (x,y,z) = (j-1, i-1, k-1)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


ras = vox2ras(1:3, 1:3) * vox + repmat(vox2ras(1:3, 4), 1, size(vox, 2));