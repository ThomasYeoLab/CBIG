function corr_mat = CBIG_self_corr(s_series)

% Calculate self correlation within a matrix.
%
%     corr_mat = CBIG_self_corr(s_series) 
%     Inputs:
%         s_series: T x N matrix, T is number of frames, N is number of voxels
%     Outputs:
%         corr_mat: N x N matrix, N is number of voxels
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


s_series = bsxfun(@minus, s_series, mean(s_series, 1));
s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

corr_mat = s_series' * s_series;
