function corr_mat = CBIG_nancorr(s_series, t_series)

% Calculate correlation matrix between each column of two matrix where
% there could be nan values in the matrices.
% 
%     corr_mat = CBIG_nancorr(s_series, t_series)
%     Input:
%         s_series: T x N1 matrix, T is number of frames, N1 number of voxels
%         t_series (optional): T x N2 matrix, T is number of frames, N2
%         number of voxels. If t_series is not passed in, then the function
%         does a self corr with s_series itself.
%     Output:
%         corr_mat: N1 x N1 matrix (self-corr) or N1 x N2 matrix (otherwise)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

counts = (sum(~isnan(s_series),1) - 1).^(0.5);
s_series = bsxfun(@minus, s_series, CBIG_nanmean(s_series, 1));
s_series = bsxfun(@times, s_series, 1./(counts.*CBIG_nanstd(s_series,0,1)));

nanind_s = isnan(s_series);

s_series(nanind_s) = 0;

if exist('t_series', 'var')
    countt = (sum(~isnan(t_series),1) - 1).^(0.5);
    t_series = bsxfun(@minus, t_series, CBIG_nanmean(t_series, 1));
    t_series = bsxfun(@times, t_series, 1./(countt.*CBIG_nanstd(t_series,0, 1)));
    
    nanind_t = isnan(t_series);
  
    t_series(nanind_t) = 0;

    corr_mat = s_series' * t_series;
else
    corr_mat = s_series' * s_series;
end
