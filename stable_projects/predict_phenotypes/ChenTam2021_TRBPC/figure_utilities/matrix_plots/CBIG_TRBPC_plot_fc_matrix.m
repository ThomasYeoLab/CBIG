function CBIG_TRBPC_plot_fc_matrix(corr_mat, scalelim, cmap, filename_prefix)

% CBIG_TRBPC_plot_fc_matrix(corr_mat, scalelim, cmap, filename_prefix)
% 
% This function simply plots a 419x419 matrix (that is already arranged into
% the major networks and subcortical regions).
% 
% Input:
%     - corr_mat: a 419x419 matrix
%     - scalelim: a vector of length 2, containing integers or floats for the 
%         minimum and maximum values for the desired color scale, e.g. [-1 1]
%     - cmap: a string for the name of the diesired color map;
%         acceptable strings are: 'parula', 'hot_cold', 'hot'
%     - filename_prefix: a string to name your output (optional, if
%         provided, a .png will be saved
%
% Written by Siyi Tang, Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% function starts here

% Load colormap
load corr_mat_colorscale.mat;

% Get grids
[major_grid, minor_grid, subcor_grid] = LabelsRearrangebyNetwork;

% Plot corr_mat using imagesc
figure
imagesc(corr_mat)

if ((nargin < 3) || (isempty(cmap)))
    cmap = 'parula';
end
if strcmp(cmap, 'hot_cold')
    set(gcf,'Colormap',rbmap2);
elseif strcmp(cmap, 'hot')
    colormap('hot')
end
    
% Generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));

xlim(gca,[1 size(corr_mat, 1)]);
ylim(gca,[1 size(corr_mat, 1)]);
postn = get(gca, 'Position');
postn(2) = 0.15;
set(gca, 'Position', postn);

% Set colorbar
hcol=colorbar('peer',gca,'SouthOutside','FontSize',15);
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/4; % Halve the thickness
cpos(3)=cpos(3)*0.75; % Reduce length
cpos(1)=cpos(1) + 0.1; % Move it to the center
cpos(2)=cpos(2) - 0.06; % Move it down outside the plot
set(hcol,'Position',cpos);


% Set color limit 
if ((nargin < 2) || (isempty(scalelim)))
	collim = max(max(abs(corr_mat)));
	scalelim = [-1*collim, 1*collim];
end

set(gca, 'CLim', scalelim);

axis equal;
grid off;
axis([-5 size(corr_mat, 1)+5.5 -5 size(corr_mat, 1)+5.5]);
set(gca, 'Visible', 'off')
set(gcf, 'color', 'white');

% Generate major and minor gridlines 
patch(xline(:,subcor_grid), yline(:,subcor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.2);
patch(yline(:,subcor_grid), xline(:,subcor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.2);
patch(xline(:,minor_grid), yline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(xline(:,major_grid), ymaj(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.1);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.1);
    
% save figure
if ((nargin == 4) && ~isempty(filename_prefix))
    filenamefin = [filename_prefix '_minsc' num2str(scalelim(1), '%.1e') '_maxsc' num2str(scalelim(2), '%.1e') '.png'];
    set(gca,'FontSize',15);
    saveas(gcf,filenamefin);
    %remove_white_board(filenamefin,filenamefin);
    % save as .eps first, and then .png
%     hgexport(gcf, filenamefin);
%     eps2xxx([filenamefin '.eps'], {'png'});

    close all
end


function [x, y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
ymaj = repmat(ymaj,1,(n-1));


function [major_grid, minor_grid, subcor_grid] = LabelsRearrangebyNetwork

% load original cortical networks labels
networks_path = fullfile(getenv('CBIG_CODE_DIR'),...
    'stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/');
lh_annot_file = [networks_path 'lh.Schaefer2018_400Parcels_17Networks_order.annot'];
rh_annot_file = [networks_path 'rh.Schaefer2018_400Parcels_17Networks_order.annot'];

[lh_vertex_labels, lh_colortable] = CBIG_read_annotation(lh_annot_file);
[rh_vertex_labels, rh_colortable] = CBIG_read_annotation(rh_annot_file);

lh_label = lh_colortable.struct_names(2:end);
rh_label = rh_colortable.struct_names(2:end);

% hard-coded, assuming the original input subcortical labels are in ascending order
subcor_labelname = {'Cerebellum_Left';'Thalamus_Left';'Caudate_Left';'Putamen_Left';'Pallidum_Left';...
    'BrainStem';'Hippocampus_Left';'Amygdala_Left';'Accumbens_Left'; ...
    'DiencephalonVentral_Left';'Cerebellum_Right';'Thalamus_Right';'Caudate_Right';'Putamen_Right';...
    'Pallidum_Right';'Hippocampus_Right';'Amygdala_Right';'Accumbens_Right';'DiencephalonVentral_Right'};

major_grid = [];
minor_grid = [];
subcor_grid = [];
%major_acc_index = [1, 4, 7, 8, 10, 12, 14, 16];
major_acc_index = [1, 4, 7, 9, 11, 13, 15, 17];
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
tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'LimbicA'; 'LimbicB';...
    'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'; ...
    'Accumbens'; 'Caudate'; 'Pallidum'; 'Putamen'; 'Thalamus'; 'Amygdala'; 'Hippocampus'; 'BrainStem';...
    'DiencephalonVentral'; 'Cerebellum'};
%tmplate2 = {'TempPole'; 'OFC'};

newlabel = [];

curInd = 0;
for j = 1:numel(tmplate)
    ind = find(strcmp(all_netw,tmplate(j)));
    
%    % for Limbic networks, further separate the networks into TempPole and OFC
%    if sum(strcmp('Limbic', tmplate(j))) ~= 0
%        ind2 = [];
%        for s = 1:numel(tmplate2)
%            if sum(strcmp(tmplate2(s), all_subnet(ind))) ~= 0
%                ind2 = [ind2; ind(find(strcmp(tmplate2(s), all_subnet(ind))))];
%                minor_grid = [minor_grid curInd+size(ind2,1)];
%            end
%            
%        end
%        if numel(ind) ~= numel(ind2)
%            disp('Wrong Index')
%        end
%        ind = ind2;
%    end

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
