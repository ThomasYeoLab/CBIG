function [corr, frac_diff] = CBIG_preproc_CensorCorrAndFracDiff(vol1, vol2)

% [corr, frac_diff] = CBIG_preproc_CensorCorrAndFracDiff(vol1, vol2)
%
% This function computes the correlation and fractional difference for each
% voxel between two N x T volumetric data matrices vol1 and vol2. 
% The fractional difference is computed by:
% 1. Sum the absolute difference between vol1 and vol2 along the time
%    dimension for each voxel.
% 2. Sum the absolute intensity of vol1 along the time dimension for each
%    voxel.
% 3. Divide the result of step 1 by the result of step 2 (element-wise
%    division).
%
% Input: 
%     - vol1:
%       a N x T matrix, where N is the number of voxels, T is the length of
%       timeseries.
%     
%     - vol2:
%       a N x T matrix, where N is the number of voxels, T is the length of
%       timeseries.
%
% Output:
%     - corr:
%       a N x 1 vecotr, entry corr(k) is the correlation between vol1(k, :)
%       and vol2(k, :).
%
%     - frac_diff:
%       a N x 1 vector, entry frac_diff(k) is the fractional difference
%       between vol1(k, :) and vol2(k, :).
%
% Example:
% fmri1 = MRIread(nifti_file1);
% fmri2 = MRIread(nifti_file2);
% vol1 = reshape(fmri1.vol, size(fmri1.vol,1)*size(fmri1.vol,2)*size(fmri1.vol,3), size(fmri1.vol,4));
% vol2 = reshape(fmri2.vol, size(fmri2.vol,1)*size(fmri2.vol,2)*size(fmri2.vol,3), size(fmri2.vol,4));
% [corr, frac_diff] = CBIG_preproc_CensorCorrAndFracDiff(vol1, vol2)
%
% Date: Jun.14, 2016
%
% Written by Jingwei Li.
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% each row of vol1/vol2 is one time series

%correlation
vol1_corr = bsxfun(@minus, vol1, mean(vol1,2)); 
vol1_corr = bsxfun(@times, vol1_corr, 1./sqrt(sum(vol1_corr.^2,2)));
vol2_corr = bsxfun(@minus, vol2, mean(vol2,2));
vol2_corr = bsxfun(@times, vol2_corr, 1./sqrt(sum(vol2_corr.^2,2)));

corr = sum(vol1_corr .* vol2_corr, 2);


%fractional difference
sum_abs = sum(abs(vol1-vol2), 2);
frac_diff = sum_abs ./ sum(abs(vol1), 2);


end