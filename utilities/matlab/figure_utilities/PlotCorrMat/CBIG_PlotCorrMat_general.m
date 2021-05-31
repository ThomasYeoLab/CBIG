function CBIG_PlotCorrMat_general(res, corr_mat, parcelname, network_order_struct, scalelim)

% CBIG_PlotCorrMat_general(res, corr_mat, parcelname, network_order_struct, scalelim)
%
% This function is a general script which visualizes ROI2ROI correlation matrix of any parcellation.
% The correlation matrix <corr_mat> can contain subcortical ROIs.
% The original ROI order of <corr_mat> is defined by <parcelname>.
% The ROIs will be visualized with the network orderring which is defined as cell structure
% <network_order_struct>:
%  {{'Net1','Net2'}, %thick lines
%   {'Net3'}, %thick lines
%   {'Net4'}}
% Thin white lines will be drawn between different subnetworks. Thick white lines will be drawn 
% between different networks based on above <network_order_struct>. For the example above, the 
% white lines will be drawn as follows:
% Net1
% ------ thin line
% Net2
% ====== thick line
% Net3
% ====== thick line
% Net4
%
% Input:
%   - res: (scalar)
%       The number of CORTICAL ROIs. For example:
%       res = 400 (not 419!) for Schaefer400+19subcortical
%       res = 100 for Schaefer100 without subcortical
%       res = 114 for Yeo2011 17 networks with split components.
%
%   - corr_mat: (KxK matrix)
%       The ROI2ROI correlation matrix. K = res if there is no subcortical ROI, else K = res+M if
%       there are M subcortical ROIs. The corr_mat is assumes to be formated as follows:
%       If there is no subcortical ROI:
%       corr_mat = [ lh2lh lh2rh;
%                    rh2lh rh2rh; ]
%       If there are subcortical ROIs:
%       corr_mat = [ lh2lh       lh2rh       lh2subcor;
%                    rh2lh       rh2rh       rh2subcor;
%                    subcor2lh   subcor2rh   subcor2subcor;]                    
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
%   - scalelim: ([lim_min limmax])
%       The range of the colormap. The <scalelim> should be a 1x2 vector where <lim_min> and <lim_max>
%       are the mininum and maximum value of the range. By default, <lim_min> is the -1*maximum
%       abosulute value of corr_mat, <lim_max> is the maximum abosulute value of corr_mat.
%
% Examples:
% corr_mat = rand(6,6)-0.5;
% parcelname = {'Net3_componnet1';
%               'Net1_componnet1';
%               'Net2_componnet1';
%               'Net1_componnet2';
%               'Net2_componnet2';
%               'Net4_componnet1';}
% network_order_struct = {{'Net1','Net2'}, %thick lines
%                         {'Net3'}, %thick lines
%                         {'Net4'}} 
% CBIG_PlotCorrMat_general(6, corr_mat, parcelname, network_order_struct, [-0.5 0.5]);
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(size(corr_mat,1) > res)
    disp(['Dimension of corr_mat is higher than res = ' num2str(res)])
    disp('Visualizing FC matrix with subcortical regions ...')
    disp('NOTE: the subcortical regions should be the last few ROIs in the input correlation matrix.')
    disp('NOTE: we will not reorder the subcortial ROIs.')
    subcortical_flag = 1;
elseif(size(corr_mat,1) == res)
    subcortical_flag = 0;
else
    error(['Dimension of corr_mat is higher than res = ' num2str(res)]);
end

%% set color limit
if ((nargin < 5) || (isempty(scalelim)))
    collim = max(max(abs(corr_mat)));
    scalelim = [-1*collim, 1*collim];
end

[Index, major_grid, minor_grid,network_order_name] = CBIG_PlotCorrMat_reorder_labels(res,...
    parcelname, network_order_struct);
% get re-indexing of ROIs
if(subcortical_flag == 1)
    Index = [Index; ((res+1):size(corr_mat,1))'];
    major_grid = [major_grid res];
end
corr_mat = corr_mat(Index,Index);

CBIG_PlotCorrMatNoOrder(corr_mat, scalelim);

%% generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));

%% generate major and minor grid lines

% cortical sub-networks
patch(xline(:,minor_grid), yline(:,minor_grid),'w', ...
'edgecolor', 'w', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', ...
'edgecolor', 'w', 'Linewidth', 0.2, 'EdgeAlpha', 0.6);

% cortical major networks
patch(xline(:,major_grid), ymaj(:,major_grid),'w', ...
'edgecolor', 'w', 'Linewidth', 1.5,'EdgeAlpha', 0.9);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', ...
'edgecolor', 'w', 'Linewidth', 1.5,'EdgeAlpha', 0.9);

% add network name
for netname = 1:length(network_order_name)
    if(netname == 1)
        ypos = xline(1,minor_grid(netname))./2;
    elseif(netname == length(network_order_name))
        ypos = xline(1,minor_grid(netname-1)) + (res + 0.5 - xline(1,minor_grid(netname-1)))./2;
    else
        ypos = xline(1,minor_grid(netname-1)) + (xline(1,minor_grid(netname)) - ...
             xline(1,minor_grid(netname-1)))./2;
    end
    text(0,ypos,network_order_name{netname},'HorizontalAlignment','right');
end
if(subcortical_flag == 1)
    ypos = res + (size(corr_mat,1) - res)./2;
    text(0,ypos,'Subcortex','HorizontalAlignment','right');
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