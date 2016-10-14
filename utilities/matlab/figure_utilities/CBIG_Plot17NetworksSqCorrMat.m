function CBIG_Plot17NetworksSqCorrMat(fig_num, corr_mat, clow, chigh, color_map)

% CBIG_Plot17NetworksSqCorrMat(fig_num, corr_mat, clow, chigh, color_map)
%
% This function is used to visualize network correlation map. Use axis 
% lines with equal lengths. Adjust the increments between data units 
% accordingly. Values less than or equal to clow map to the first color in 
% the color_map and values greater than or equal to chigh map to the last 
% color in the color_map.
%
% Input:
%      -fig_num:
%       Output figure number
%
%      -corr_mat:
%       NxM (N can euqal to M) ROIs2ROIs correlation matrix, N is the
%       number of ROIs in xhemi, M is the number of ROIs in yhemi.
%
%
%      -clow, chigh:
%       Values less than or equal to clow will map to the first color in 
%       the color_map and values greater than or equal to chigh will map to 
%       the last color in the color_map. Defaut clow is -1, default chigh
%       is 1.
%
%      -color_map:
%       Kx3 color matrix. K is the number of different colors. Each row of
%       this matrix corresponds to the rgb value. There is a default
%       colormap from dark black to light red.
%
% Example:
% CBIG_Plot17NetworksSqCorrMat(1, correlation_mat)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(~exist('color_map', 'var'))
    cyan_black = [linspace(0, 0, 32)' linspace(255, 0, 32)' linspace(255, 0, 32)'];
    black_red   = [linspace(0, 255, 32)' linspace(0, 0, 32)' linspace(0, 0, 32)'];
    color_map = [cyan_black; black_red]/255;
end

if(~exist('clow', 'var'))
    clow = -1;
end

if(~exist('chigh', 'var'))
    chigh = 1;
end

h = figure(fig_num); gpos = get(h, 'Position');
gpos(3) = 1150; gpos(4) = 1150; set(h, 'Position', gpos);
imagesc(corr_mat, [clow chigh]);
set(gca, 'YDir', 'normal');
colormap(color_map); 
axis square; axis off;
