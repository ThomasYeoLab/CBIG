function CBIG_PlotCorrMatNoOrder(corr_mat, scalelim)

% CBIG_PlotCorrMatNoOrder(corr_mat, scalelim)
%
% This script visualize the ROI2ROI correlation matrix <corr_mat> with given colormap range <scalelim>.
% The script will not perform and reordering but simply visualize correlation matrix as the original
% order.
%
%   - corr_mat: (KxK matrix)
%       The ROI2ROI correlation matrix. 
%
%   - scalelim: ([lim_min limmax])
%       The range of the colormap. The <scalelim> should be a 1x2 vector where <lim_min> and <lim_max>
%       are the mininum and maximum value of the range. By default, <lim_min> is the -1*maximum
%       abosulute value of corr_mat, <lim_max> is the maximum abosulute value of corr_mat.
%
% Examples:
% corr_mat = rand(6,6)-0.5;
% CBIG_PlotCorrMatNoOrder(corr_mat, [-0.5 0.5]);  
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% load colormap
load('corr_mat_colorscale.mat');

%% plot correlation matrix using imagesc
figure; imagesc(corr_mat); set(gcf,'Colormap',rbmap2);

xlim(gca,[1 size(corr_mat, 1)]);
ylim(gca,[1 size(corr_mat, 1)]);
postn = get(gca, 'Position');
postn(2) = 0.15;
set(gca, 'Position', postn);

%% set colorbar
hcol=colorbar('peer',gca,'SouthOutside','FontSize',15);
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/4; % Halve the thickness
cpos(3)=cpos(3)*0.75; % Reduce length
cpos(1)=cpos(1) + 0.1; % Move it to the center
cpos(2)=cpos(2) - 0.12; % Move it down outside the plot
set(hcol,'Position',cpos);

%% set color limit
if ((nargin < 2) || (isempty(scalelim)))
    collim = max(max(abs(corr_mat)));
    scalelim = [-1*collim, 1*collim];
end

set(gca, 'CLim', scalelim);
if(~strcmp(version(), '8.3.0.532 (R2014a)'))
    hcol.Ticks = linspace(scalelim(1), scalelim(2), 7);
    hcol.TickLabels = num2cell(hcol.Ticks);
    hcol.FontSize = 8;
else
    set(hcol, 'XTick', linspace(scalelim(1), scalelim(2), 7));
    set(hcol, 'XTickLabel', num2str(linspace(scalelim(1), scalelim(2), 7)'));
    set(hcol, 'FontSize', 8);
end

axis equal;
grid off;
axis([-5 size(corr_mat, 1)+5.5 -5 size(corr_mat, 1)+5.5]);
set(gca, 'Visible', 'off')
set(gcf, 'color', 'white');
set(gcf, 'Position', [1 1 700 700]);