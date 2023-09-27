function corr_mat = CBIG_corr(s_series, t_series)

% Calculate correlation matrix between each column of two matrix.
% 
%     corr_mat = CBIG_corr(s_series, t_series)
%     Input:
%         s_series: T x N1 matrix, T is number of frames, N1 number of voxels
%         t_series: T x N2 matrix, T is number of frames, N2 number of voxels
%     Output:
%         corr_mat: N1 x N2 matrix
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


s_series = bsxfun(@minus, s_series, mean(s_series, 1));
s_series = bsxfun(@times, s_series, 1./sqrt(sum(s_series.^2, 1)));

t_series = bsxfun(@minus, t_series, mean(t_series, 1));
t_series = bsxfun(@times, t_series, 1./sqrt(sum(t_series.^2, 1)));

corr_mat = s_series' * t_series;
