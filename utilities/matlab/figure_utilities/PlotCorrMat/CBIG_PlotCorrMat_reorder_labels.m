function [Index, major_grid, minor_grid, network_order_name] = CBIG_PlotCorrMat_reorder_labels(res,...
    parcelname, network_order_struct)

% [Index, major_grid, minor_grid, network_order_name] = CBIG_PlotCorrMat_reorder_labels(res,...
% parcelname, network_order_struct)
%
% This script reads in the original ROI order and the visualization ROI order, it will find the location
% to draw thick and thin white grids to separate network structure. The subcortical part will be separated
% from cortical part by a thick white line. This script will not consider any subcortical ROI name.
%
% Input:
%   - res: (scalar)
%       The number of CORTICAL ROIs. For example:
%       res = 400 (not 419!) for Schaefer400+19subcortical
%       res = 100 for Schaefer100 without subcortical
%       res = 114 for Yeo2011 17 networks with split components.                
%
%   - parcelname: (resx1 cell vector)
%       The original order of cortical ROIs. For example:
%       parcelname = {'Net3_componnet1';
%                     'Net1_componnet1';
%                     'Net2_componnet1';
%                     'Net1_componnet2';
%                     'Net2_componnet2';
%                     'Net4_componnet1';}
%
%   - network_order_struct: (cell vector contains sub-cell elements)
%       The visualization order of cortical ROIs. Thick white grids will be drawn between sub-cells.
%       For example:
%       {{'Net1','Net2'}, %thick lines
%        {'Net3'}, %thick lines
%        {'Net4'}}
%
% Output:
%   - Index: (resx1 vector)
%       The reordering indices.
%
%   - major_grid: (Kx1 vector)
%       The locations to draw thick white grids. For example, if major_grid = [4,5],
%       the thick white grids will be drawn after the 4th and 5th reordered ROI.
%       The thick white grids are used to separate the major network structure.
%
%   - minor_grid: (Mx1 vector)
%       The locations to draw thin white grids. For example, if major_grid = [2,4,5],
%       the thin white grids will be drawn after the 2th, 4th and 5th reordered ROI.
%       The thin white grids are used to separate different networks.
%
% Examples:
% parcelname = {'Net3_componnet1';
%               'Net1_componnet1';
%               'Net2_componnet1';
%               'Net1_componnet2';
%               'Net2_componnet2';
%               'Net4_componnet1';}
% network_order_struct = {{'Net1','Net2'}, %thick lines
%                         {'Net3'}, %thick lines
%                         {'Net4'}}  
% [Index, major_grid, minor_grid, network_order_name] = CBIG_PlotCorrMat_reorder_labels(6,...
%    parcelname, network_order_struct);
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

net_index = 0;
major_acc_index = [];
for no = 1:length(network_order_struct)
    net_index = net_index + length(network_order_struct{no});
    major_acc_index = [major_acc_index net_index];

    if (no == 1)
        network_order_name = network_order_struct{no}; 
    else
        network_order_name = [network_order_name network_order_struct{no}]; 
    end
end
% cortical
seg = [];
for i = 1:res
    netw = split(parcelname{i}, '_');
    if(i == 1)
        for t = 1:length(network_order_name)
            check_seg = ismember(netw, network_order_name{t});
            if(sum(check_seg) ~= 0)
                seg = find(check_seg == 1);
            end
        end
    end
    if(isempty(seg))
        error('Cannot find network_order_struct names in parcelname.');
    end
    tmp = netw(seg);
    all_netw(i,1) = tmp;
end
major_grid = [];
minor_grid = [];
newlabel = [];

curInd = 0;
for j = 1:numel(network_order_name)
    ind = find(strcmp(all_netw,network_order_name(j)));
    curInd = curInd+size(ind,1);
    if (j~=numel(network_order_name))
        minor_grid = [minor_grid curInd];
        if (any(j==major_acc_index))
            major_grid = [major_grid curInd];
        end 
    end
    newlabel = cat(1, newlabel,parcelname(ind));
end

% create indexing from old labeling to new labeling
Index = zeros(res,1);

for k = 1:res
    Index(k) = find(strcmp(newlabel(k),parcelname));
end

fin_label = parcelname(Index);
    