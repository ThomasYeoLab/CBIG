function CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(lh2lh_corrmat, ...
    lh2rh_corrmat, rh2rh_corrmat, lh2subcor_corrmat, rh2subcor_corrmat, ...
    subcor2subcor_corrmat, scalelim, filename_prefix)

% CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(lh2lh_corrmat, ...
%     lh2rh_corrmat, rh2rh_corrmat, lh2subcor_corrmat, rh2subcor_corrmat, ...
%     subcor2subcor_corrmat, scalelim, filename_prefix)
%
% This function draws the 419x419 correlation matrix (400 cortical ROIs +
% 19 subcortical ROIs). 
% Major networks are separated by thick white grid lines.
% Thin white grid lines separate the breakdowns of major networks and the
% 19 subcortical regions.
% Subcortical structures in the striatum are arranged together in the plot.
% Details on the ordering of the cortical networks and subcortical structures
% are shown below.
%
% This function assumes that the order of the subcortical ROIs in the
% input correlation matrices is in ascending order based on the labels in
% $FREESURFER_HOME/ASegStatsLUT.txt file.
%
% Ordering of major networks and subcortical regions from left to right:
%
%       Major network/Subcortical regions        Sub-network
%     
%     Cortical networks:
%
%     1)  Default:                                TempPar
%                                                 DefaultC 
%                                                 DefaultB
%                                                 DefaultA
%
%     2)  Control:                                ContC
%                                                 ContB 
%                                                 ContA
%
%     3)  Limbic
%
%     4)  SalVentAttn:                            SalVentAttnB
%                                                 SalVentAttnA
%
%     5)  DorsAttn:                               DorsAttnB
%                                                 DorsAttnA
%
%     6)  SomMot:                                 SomMotB
%                                                 SomMotA
%
%     7)  Visual:                                 VisPeri
%                                                 VisCent
%     
%     Subcortical structures:
%     The following 4 striatum structures are arranged together:
%     1)  Accumbens
%     2)  Caudate
%     3)  Pallidum
%     4)  Putamen
%
%     5)  Thalamus
%
%     6)  Amygdala
%     
%     7)  Hippocampus
%
%     8)  Brain Stem
%
%     9)  Diencephalon (Ventral)
%
%    10)  Cerebellum 
%
% Within each sub-network or subcortical region, correlation matrix entries
% start from left hemisphere, then right hemisphere entries (from left to
% right).
%
% Note that the highly "dense" white lines will become more appropriate
% when the figure is saved.
% 
% Input:
%     - lh2lh_corrmat:
%       200x200 left hemisphere cortical to left hemisphere cortical
%       regions correlation matrix
%     - lh2rh_corrmat:
%       200x200 left hemisphere cortical to right hemisphere cortical
%       regions correlation matrix
%     - rh2rh_corrmat:
%       200x200 right hemisphere cortical to right hemisphere cortical
%       regions correlation matrix
%     - lh2subcor_corrmat:
%       200x19 left hemisphere cortical to subcortical regions correlation
%       matrix. 
%       The order of the subcortical ROIs in the input correlation
%       matrix is assumed to be in ascending order based on the labels in
%       $FREESURFER_HOME/ASegStatsLUT.txt file.
%     - rh2subcor_corrmat:
%       200x19 right hemisphere cortical to subcortical regions correlation
%       matrix
%     - subcor2subcor_corrmat:
%       19x19 subcortical to subcortical regions correlation matrix
%     - scalelim:
%       Min and max scale limit 
%       If not specified, or specified as [], scale limit is 
%       -1*max(abs(corr_mat)) to max(abs(corr_mat)).
%     - filename_prefix:
%       Prefix of the saved file name 
%       If not specified, figure will not be saved.
%
% Example:
% CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(lh2lh, ...
%     lh2rh, rh2rh, lh2subcor, rh2subcor, subcor2subcor, [], 'corrmat') 
% The max/min scales depend on the maximum absolute value of the correlation matrix. 
% Save figure as 'corrmat_minsc-1_maxsc1.jpg'.
%
% CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid(lh2lh, ...
%     lh2rh, rh2rh, lh2subcor, rh2subcor, subcor2subcor, [-0.5 0.5])
% Will not save figure.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% rearrange input correlation matrices into a 419x419 correlation matrix
corr_mat = [lh2lh_corrmat lh2rh_corrmat lh2subcor_corrmat; ...
    lh2rh_corrmat' rh2rh_corrmat rh2subcor_corrmat; ...
    lh2subcor_corrmat' rh2subcor_corrmat' subcor2subcor_corrmat];
corr_mat = (corr_mat+corr_mat')/2;
[Index, major_grid, minor_grid, subcor_grid] = CBIG_LabelsRearrangebyNetwork;
% get re-indexing of ROIs
corr_mat = corr_mat(Index,Index);

%% load colormap
load('corr_mat_colorscale.mat');

%% plot correlation matrix using imagesc
figure; imagesc(corr_mat); set(gcf,'Colormap',rbmap2);

%% generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));

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
if ((nargin < 7) || (isempty(scalelim)))
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

%% generate major and minor grid lines
% subcortical
patch(xline(:,subcor_grid), yline(:,subcor_grid),'w', ...
    'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.2);
patch(yline(:,subcor_grid), xline(:,subcor_grid),'w', ...
    'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.2);

% cortical sub-networks
patch(xline(:,minor_grid), yline(:,minor_grid),'w', ...
    'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', ...
    'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);

% cortical major networks
patch(xline(:,major_grid), ymaj(:,major_grid),'w', ...
    'edgecolor', 'w', 'Linewidth', 1.1);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', ...
    'edgecolor', 'w', 'Linewidth', 1.1);

%% save figure
if ((nargin == 8) && ~isempty(filename_prefix))
    filenamefin = [filename_prefix '_minsc' num2str(scalelim(1), '%.2f') '_maxsc' num2str(scalelim(2), '%.2f') '.jpg'];
    set(gcf, 'PaperPositionMode', 'auto','Renderer', 'ZBuffer');
    tmp = get(gcf, 'PaperPosition');
    set(gcf, 'PaperPosition', [0 0 tmp(3)*3563/3500 tmp(4)*2581/2625]); 
    print(gcf, '-djpeg', '-r600', filenamefin);
    close all
end

end

%% sub-function to generate grid lines
function [x, y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
% For short grid (same as the figure boundary) ymaj = [ 0.5 n+0.5 nan]';
ymaj = repmat(ymaj,1,(n-1));

end

