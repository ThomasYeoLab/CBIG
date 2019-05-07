function CBIG_Plot_Schaefer400_17NetworksRearrCorrMat_WhiteGrid(lh2lhcorrmat, ...
    lh2rhcorrmat, rh2rhcorrmat, scalelim, filename_prefix)

% CBIG_Plot_Schaefer400_17NetworksRearrCorrMat_WhiteGrid(lh2lhcorrmat, ...
%    lh2rhcorrmat, rh2rhcorrmat, scalelim, filename_prefix)
% Example : CBIG_Plot17NetworksRearrCorrMat_WhiteGrid(lhlh, lhrh, rhrh, [], 'corrmat') 
%               Max scale depends on maximum of correlation matrix value
%               Save figure as corrmat_minsc-1_maxsc1.jpg
%
%           CBIG_Plot17NetworksRearrCorrMat_WhiteGrid(lhlh, lhrh, rhrh, [-0.5 0.5])
%               Will not save figure
%
% Note that the highly "dense" white lines will become more appropriate when the figure is saved.
%
% inputs are 	1) LeftHemisphere to LeftHemisphere correlation matrix
%               2) LeftHemisphere to RightHemisphere correlation matrix
%               3) RightHemisphere to RightHemisphere correlation matrix
%               4) Scale limit : min and max scale limit
%                  If not specified, or specified as [], 
%                  scale limit will be -1*max(abs(corr_mat)) to max(abs(corr_mat))
%               5) Filename : basename. 
%                  E.g. 'RWRestingCorrMat', final filename = RWRestingCorrMat_min-3_max3.jpg
%                  If not specified, figure will not be saved
%
% Major networks separated by thick white grid lines
% Thin white grid lines separate the breakdowns of major networks
%
% Ordering of major networks from left to right 
%
%       Major network     Sub-network
%     
%     1)  Default :       TempPar
%                         DefaultC
%                         DefaultB
%                         DefaultA
%         
%     2)  Control :       ContC
%                         ContB
%                         ContA
%         
%     3)  Limbic
%     
%     4)  SalVentAttn :   SalVentAttnB
%                         SalVentAttnA
%         
%     5)  DorsAttn :      DorsAttnB
%                         DorsAttnA
%         
%     6)  SomMot :        SomMotB
%                         SomMotA
%     
%     7)  Visual :        VisPeri
%                         VisCent
%
% Within each sub-network, correlation matrix entries start from left
% hemisphere, then right hemisphere entries (from left to right).
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% Rearrange corrmat data
corr_mat = [lh2lhcorrmat lh2rhcorrmat; lh2rhcorrmat' rh2rhcorrmat];
corr_mat = (corr_mat+corr_mat')/2;
% Index = open('indexList.mat');
% Index = Index.index;
[Index, major_grid, minor_grid] = LabelsRearrangebyNetwork;
% % Get re-indexing of ROIs
% save('idxnew.mat', 'Index')
corr_mat = corr_mat(Index,Index);


% Load colormap
load corr_mat_colorscale.mat;


% Plot corr_mat using imagesc
figure; imagesc(corr_mat); set(gcf,'Colormap',rbmap2);
    
% Generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));
%patch(xline, yline,'w', 'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.1)
%patch(yline, xline,'w', 'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.1)

xlim(gca,[1 size(corr_mat, 1)]);
ylim(gca,[1 size(corr_mat, 1)]);
postn = get(gca, 'Position');
postn(2) = 0.15;
set(gca, 'Position', postn);

% Set colorbar
hcol=colorbar('peer',gca,'SouthOutside');
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/4; % Halve the thickness
cpos(3)=cpos(3)*0.75; % Reduce length
cpos(1)=cpos(1) + 0.1; % Move it to the center
cpos(2)=cpos(2) - 0.12; % Move it down outside the plot
set(hcol,'Position',cpos);

% Set color limit 
if ((nargin < 4) || (isempty(scalelim)))
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

% Generate major and minor gridlines 
  
% 	major_grid = [16 95 156 180 231 283 353];
%   minor_grid = [16 29 61 95 107 132 156 169 180 197 231 256 283 314 353 376];
patch(xline(:,minor_grid), yline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(xline(:,major_grid), ymaj(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.1);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.1);
    
% save figure
if ((nargin == 5) && ~isempty(filename_prefix))
    filenamefin = [filename_prefix '_minsc' num2str(scalelim(1), '%.2f') '_maxsc' num2str(scalelim(2), '%.2f') '.jpg'];
    set(gcf, 'PaperPositionMode', 'auto','Renderer', 'ZBuffer'); 
    print(gcf, '-djpeg', '-r600', filenamefin); 
    % Crop image
%     origfile = imread(fullfile(filenamefin));
%     croppedimg = imcrop(origfile, [750 50 2350 2700]);
%     imwrite(croppedimg, fullfile(filenamefin));
    close all
end


function [x, y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
% For short grid (same as the figure boundary)
% ymaj = [ 0.5 n+0.5 nan]';
ymaj = repmat(ymaj,1,(n-1));



function [Index, major_grid, minor_grid] = LabelsRearrangebyNetwork

% Load original labels
networks_path = fullfile(getenv('CBIG_CODE_DIR'), ...
    'stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/');
lh_annot_file = [networks_path 'lh.Schaefer2018_400Parcels_17Networks_order.annot'];
rh_annot_file = [networks_path 'rh.Schaefer2018_400Parcels_17Networks_order.annot'];
[lh_vertex_labels, lh_colortable] = CBIG_read_annotation(lh_annot_file);
[rh_vertex_labels, rh_colortable] = CBIG_read_annotation(rh_annot_file);

lh_label = lh_colortable.struct_names(2:end);
rh_label = rh_colortable.struct_names(2:end);

major_grid = [];
minor_grid = [];
major_acc_index = [1, 4, 7, 8, 10, 12, 14];

%% 
% jwpath = fullfile('share', 'users', 'imganalysis', 'yeolab', 'data', ...
%     'PublicParcellations', 'Schaefer', '400_parcels', 'fsaverage6', ...
%     'annot', 'lh.Schaefer2016_400Parcels_17Networks_colors_23_05_16.annot');
% % [jwlabel, label, colorlabel] = read_annotation(jwpath)
%  
% lh_label = colorlabel.struct_names(2:201);
% rh_label = colorlabel.struct_names(202:401);

%% 



% labeldir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference');
% 
% lh_fid = fopen(fullfile(labeldir, 'lh.Yeo2011_17Networks_N1000.split_components.txt'));
% rh_fid = fopen(fullfile(labeldir, 'rh.Yeo2011_17Networks_N1000.split_components.txt'));
% 
% lh_label = textscan(lh_fid, '%s %s %s', 'delimiter', '.'); lh_label = lh_label{1,2}; 
% rh_label = textscan(rh_fid, '%s %s %s', 'delimiter', '.'); rh_label = rh_label{1,2}; 

lhrh_label = {lh_label{:} rh_label{:}}';
lhrh_netw = cell(numel(lhrh_label),1);
lhrh_subnet = cell(numel(lhrh_label),1);

for i = 1:numel(lhrh_label);
    netw = textscan(char(lhrh_label(i,1)), '%s %s %s %s', 'delimiter', '_');
    tmp = netw{1,3};
    lhrh_netw(i,1) = tmp;
    subnet = netw{1,4};
    if ~isempty(subnet)
        lhrh_subnet(i,1) = subnet;
    end
end


% Arrange new label based on template order
tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; ...
    'ContA'; 'Limbic'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; ...
    'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'};
tmplate2 = {'TempPole'; 'OFC'};
% initiate new label
newlabel = [];
%% 

curInd = 0;
for j = 1:numel(tmplate);
    ind = strmatch(tmplate(j), lhrh_netw);
    
    % For Limbic, further separate the networks to TempPole and OFC
    if ~isempty(strmatch('Limbic', tmplate(j)));
        ind2 = [];
        for s = 1:numel(tmplate2);
            if ~isempty(strmatch(tmplate2(s), lhrh_subnet(ind)))
                ind2 = [ind2; ind(strmatch(tmplate2(s), lhrh_subnet(ind)))];
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
        minor_grid = [minor_grid curInd];
    end
    if (any(j==major_acc_index))
        major_grid = [major_grid curInd];
    end
    newlabel = cat(1, newlabel,lhrh_label(ind));
end
%% 


% Create indexing from old labeling to new labeling
Index = zeros(numel(lhrh_label),1);

for k = 1:numel(lhrh_label);
Index(k) = strmatch(newlabel(k), lhrh_label, 'exact');
end

fin_label = lhrh_label(Index);
% Save new label
fid = fopen('networks_label_sorted_NEWORDERING.txt', 'w');
for l = 1:numel(lhrh_label);
   fprintf(fid, '%s\n', char(fin_label(l)));
end