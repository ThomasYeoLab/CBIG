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
[Index, major_grid, minor_grid, subcor_grid] = LabelsRearrangebyNetwork;
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
    print(gcf, '-djpeg', '-r600', filenamefin);
    close all
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


%% sub-function to rearrange cortical networks and subcortical regions
function [Index, major_grid, minor_grid, subcor_grid] = LabelsRearrangebyNetwork

% load original cortical networks labels
networks_path = fullfile(getenv('CBIG_CODE_DIR'), ...
    'stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/');
lh_annot_file = [networks_path 'lh.Schaefer2018_400Parcels_17Networks_order.annot'];
rh_annot_file = [networks_path 'rh.Schaefer2018_400Parcels_17Networks_order.annot'];
[lh_vertex_labels, lh_colortable] = CBIG_read_annotation(lh_annot_file);
[rh_vertex_labels, rh_colortable] = CBIG_read_annotation(rh_annot_file);

lh_label = lh_colortable.struct_names(2:end);
rh_label = rh_colortable.struct_names(2:end);

% hard-coded, assuming the original input subcortical labels are in ascending order
subcor_labelname = {'Cerebellum_Left';'Thalamus_Left';'Caudate_Left';'Putamen_Left';...
    'Pallidum_Left';'BrainStem';'Hippocampus_Left';'Amygdala_Left';'Accumbens_Left'; ...
    'DiencephalonVentral_Left';'Cerebellum_Right';'Thalamus_Right';'Caudate_Right';...
    'Putamen_Right';'Pallidum_Right';'Hippocampus_Right';'Amygdala_Right';...
    'Accumbens_Right';'DiencephalonVentral_Right'};

major_grid = [];
minor_grid = [];
subcor_grid = [];
major_acc_index = [1, 4, 7, 8, 10, 12, 14, 16];
subcor_acc_index = 17:25;

lhrh_label = {lh_label{:} rh_label{:}}';
all_label = [lhrh_label; subcor_labelname];
all_netw = cell(numel(all_label),1);
all_subnet = cell(numel(all_label),1);

for i = 1:numel(all_label)
    % cortical
    if i <= numel(lhrh_label)
        netw = textscan(char(all_label(i,1)), '%s %s %s %s', 'delimiter', '_');
        tmp = netw{1,3};
        all_netw(i,1) = tmp;
        subnet = netw{1,4};
        if ~isempty(subnet)
            all_subnet(i,1) = subnet;
        end
        
    % subcortical    
    else
        netw = textscan(char(all_label(i,1)),'%s %s', 'delimiter', '_');
        all_netw(i,1) = netw{1,1};       
    end
end

% arrange new labels based on template order
tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; ...
    'ContA'; 'Limbic'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; ...
    'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'; ...
    'Accumbens'; 'Caudate'; 'Pallidum'; 'Putamen'; 'Thalamus'; ...
    'Amygdala'; 'Hippocampus'; 'BrainStem'; 'DiencephalonVentral'; 'Cerebellum'};
tmplate2 = {'TempPole'; 'OFC'};

newlabel = [];

curInd = 0;
for j = 1:numel(tmplate)
    ind = find(strcmp(all_netw,tmplate(j)));
    
    % for Limbic networks, further separate the networks into TempPole and OFC
    if sum(strcmp('Limbic', tmplate(j))) ~= 0
        ind2 = [];
        for s = 1:numel(tmplate2)
            if sum(strcmp(tmplate2(s), all_subnet(ind))) ~= 0
                ind2 = [ind2; ind(find(strcmp(tmplate2(s), all_subnet(ind))))];
                minor_grid = [minor_grid curInd+size(ind2,1)];
            end
            
        end
        if numel(ind) ~= numel(ind2)
            disp('Wrong Index')
        end
        ind = ind2;
    end
    curInd = curInd+size(ind,1);
    if (j~=numel(tmplate))
        if (any(j==subcor_acc_index))
            subcor_grid = [subcor_grid curInd];
        else
            minor_grid = [minor_grid curInd];
        end
    end
    if (any(j==major_acc_index))
        major_grid = [major_grid curInd];
    end
    
    newlabel = cat(1, newlabel,all_label(ind));
end

% create indexing from old labeling to new labeling
Index = zeros(numel(all_label),1);

for k = 1:numel(all_label)
    Index(k) = find(strcmp(newlabel(k),all_label));
end

fin_label = all_label(Index);

% save new labels
fid = fopen('networks_label_sorted_NEWORDERING.txt', 'w');
for l = 1:numel(all_label);
    fprintf(fid, '%s\n', char(fin_label(l)));
end
