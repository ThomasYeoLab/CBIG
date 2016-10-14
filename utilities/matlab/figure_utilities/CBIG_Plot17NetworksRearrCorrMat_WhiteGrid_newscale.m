function curr_net = CBIG_Plot17NetworksRearrCorrMat_WhiteGrid_newscale(lh2lhcorrmat, lh2rhcorrmat, rh2rhcorrmat, scalelim, filename_prefix)

% curr_net = CBIG_Plot17NetworksRearrCorrMat_WhiteGrid_newscale(lh2lhcorrmat, lh2rhcorrmat, rh2rhcorrmat, scalelim, filename_prefix)
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
%                  If not specified, or specified as [], scale limit is -1*max(abs(corr_mat)) to max(abs(corr_mat))
%               5) Filename : basename. E.g. 'RWRestingCorrMat', final filename = RWRestingCorrMat_min-3_max3.jpg
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
%			  DefaultC
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
[curr_net, Index] = LabelsRearrangebyNetwork;
% % Get re-indexing of ROIs
% save('idxnew.mat', 'Index')
corr_mat = corr_mat(Index,Index);


% Load colormap
load corr_mat_colorscale_bluemod.mat;


% Plot corr_mat using imagesc
figure; imagesc(corr_mat); set(gcf,'Colormap',rbmap2);
    
% Generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));
patch(xline, yline,'w', 'edgecolor', 'w', 'Linewidth', 0.15, 'EdgeAlpha', 0.8)
patch(yline, xline,'w', 'edgecolor', 'w', 'Linewidth', 0.15, 'EdgeAlpha', 0.8)

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

axis equal;
grid off;
axis([-5 size(corr_mat, 1)+5.5 -5 size(corr_mat, 1)+5.5]);
set(gca, 'Visible', 'off')
set(gcf, 'color', 'white');

% Generate major and minor gridlines 
	minor_grid = [2 8 17 30 41 54 69 88 102 110];
	major_grid = [26 52 56 80 94 104];

	patch(xline(:,minor_grid), yline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.6, 'EdgeAlpha', 0.9)
	patch(yline(:,minor_grid), xline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.6, 'EdgeAlpha', 0.9)
	patch(xline(:,major_grid), ymaj(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.2)
	patch(ymaj(:,major_grid), xline(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.2)

% save figure
if ((nargin == 5) && ~isempty(filename_prefix))
    filenamefin = [filename_prefix '_minsc' num2str(scalelim(1), '%.2f') '_maxsc' num2str(scalelim(2), '%.2f') '.jpg'];
    set(gcf, 'PaperPositionMode', 'auto'); 
    print(gcf, '-djpeg', '-r600', filenamefin); 
    % Crop image
    origfile = imread(fullfile(filenamefin));
    croppedimg = imcrop(origfile, [750 50 2350 2700]);
    imwrite(croppedimg, fullfile(filenamefin));
    close all
end


function [x,y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -5 n+5.5 nan]';
% For short grid (same as the figure boundary)
% ymaj = [ 0.5 n+0.5 nan]';
ymaj = repmat(ymaj,1,(n-1));



function [curr_net, Index] = LabelsRearrangebyNetwork

% Load original labels
labeldir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference');

lh_fid = fopen(fullfile(labeldir, 'lh.Yeo2011_17Networks_N1000.split_components.txt'));
rh_fid = fopen(fullfile(labeldir, 'rh.Yeo2011_17Networks_N1000.split_components.txt'));
lh_label = textscan(lh_fid, '%s %s %s', 'delimiter', '.'); lh_label = lh_label{1,2}; 
rh_label = textscan(rh_fid, '%s %s %s', 'delimiter', '.'); rh_label = rh_label{1,2}; 

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
tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'Limbic'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'};
tmplate2 = {'TempPole'; 'OFC'};
% initiate new label
newlabel = [];
curr_net = NaN(size(lhrh_netw));

for j = 1:numel(tmplate);
    ind = strmatch(tmplate(j), lhrh_netw);
    curr_net(ind) = j;
    % For Limbic, further separate the networks to TempPole and OFC
    if ~isempty(strmatch('Limbic', tmplate(j)));
        ind2 = [];
        for s = 1:numel(tmplate2);
            if ~isempty(strmatch(tmplate2(s), lhrh_subnet(ind)))
                ind2 = [ind2; ind(strmatch(tmplate2(s), lhrh_subnet(ind)))];
            end
            
        end
        if numel(ind) ~= numel(ind2)
            disp('Wrong Index')
        end
        ind = ind2;
    end
    
    newlabel = cat(1, newlabel,lhrh_label(ind));
end

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

