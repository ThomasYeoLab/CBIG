function CBIG_Plot17NetworksCorrMat(fig_num, corr_mat, xhemi, yhemi, clow, chigh, color_map)

% CBIG_Plot17NetworksCorrMat(fig_num, corr_mat, xhemi, yhemi, clow, chigh, color_map)
%
% This function is used to visualize 17 network correlation map. Values 
% less than or equal to clow map to the first color in the color_map and 
% values greater than or equal to chigh map to the last color in the 
% color_map. If corr_mat is the correlation between lh and lh, then 
% xhemi = 'lh', yhemi = 'rh'.
% 
% Input:
%      -fig_num:
%       Output figure number
%
%      -corr_mat:
%       NxM (N can euqal to M) ROIs2ROIs correlation matrix, N is the
%       number of ROIs in xhemi, M is the number of ROIs in yhemi.
%
%      -xhemi, yhemi:
%       'lh' or 'rh'. xhemi and yhemi can be the same.
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
% CBIG_Plot17NetworksCorrMat(1,lhrh_corr_mat,'lh','rh', 0.1, 0.8)
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

if(~exist('xhemi', 'var'))
    xhemi = 'lh';
end

if(~exist('yhemi', 'var'))
    yhemi = 'lh';
end

xlabel_txt = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference', [xhemi '.Yeo2011_17Networks_N1000.split_components.txt']);
ylabel_txt = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference', [yhemi '.Yeo2011_17Networks_N1000.split_components.txt']);

Xlabels = textread(xlabel_txt, '%s');
for i = 1:length(Xlabels)
   label = basename(Xlabels{i});
   Xlabels{i} = label(15:end-6);
   Xlabels{i} = strrep(Xlabels{i}, '_', '\_');
end

Ylabels = textread(ylabel_txt, '%s');
for i = 1:length(Ylabels)
   label = basename(Ylabels{i});
   Ylabels{i} = label(15:end-6);
end

h = figure(fig_num); gpos = get(h, 'Position');
gpos(3) = 1150; gpos(4) = 1300; set(h, 'Position', gpos);
disp(['max val: ' num2str(max(corr_mat(:))) ', min val: ' num2str(min(corr_mat(:)))]);
imagesc(corr_mat, [clow chigh]);

set(gca, 'Position', [0.050    0.030    0.7750    0.7750]);
set(gca, 'YTick', 1:length(Ylabels));
set(gca, 'YTickLabel', Ylabels);
set(gca, 'YAxisLocation', 'right');
set(gca, 'YDir', 'normal');
set(gca, 'XTick', []);
text(1:length(Xlabels), zeros(1, length(Xlabels))+length(Ylabels)+0.75, Xlabels, 'HorizontalAlignment','left', 'rotation', 90, 'FontSize', 10, 'FontName', 'San Serif')
set(gca, 'FontSize', 10);
set(gca, 'FontName', 'San Serif');
colormap(color_map); 
colorbar('Location', 'SouthOutside');


