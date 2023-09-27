function CBIG_GraphVis_17network(data, netw_used, hemi, fileout, seed_in)

% CBIG_GraphVis_17network(data, netw_used, hemi, fileout, seed_in)
% Example usage: 
%
%   CBIG_GraphVis_17network(corr_mat, [1 3 4 6 7 9 10 11 12 15 16], 'lh', 'test', 20);
%       will plot graph with left hemisphere ROIs of networks 1 3 4 6 7 9 10 11 12 15 16
%       using seed 20, the 2D data and seed will be saved in rawdata/test.mat 
%       the figure will be saved as figures/test.eps
%   fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Yeo2011_fcMRI_clustering', '1000subjects_reference');
%
% Input :   data      -  114 x 114 matrix from 17-network parcellation with ROIs arranged based on 
%                        $CBIG_CODE_DIR/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/?h.Yeo2011_17Networks_N1000.split_components.txt
%
%           netw_used -  specifies what networks to plot. The code will print out the network you are saving, so you
%                        can experiment with what networks you want to save. Note that Limbic_TempPole and Limbic_OFC are
%                        combined as limbic, so the maximum number of networks to specify is 16. 
%                        To plot everything, netw_used = 1:16. By default netw_used = [1 3 4 6 7 9 10 11 12 15 16]; 
%
%           hemi      -  'lh', 'rh', or 'lhrh'. By default = 'lh', which plots only left hemisphere
%
%           fileout   -  prefix for the output files. 2 folders will be generated : 'figures' and 'rawdata'
%
%           seed_in   -  optional seed for random number generator, if you wish to regenerate previously
%                        generated graph, pass in the seed from that graph (from the rawdata/*.mat file)
%                        You might want to try different seeds and settle for the best results
%                        Note that at least for the 1000 subjects 114 x 114 matrix, the results are the same regardless of
%                        seed (save for a rotation factor). For example, in
%                        my data, seed 20 generates same results as seed 5 save for
%                        a rotation factor, but seed 20 results look visually more compelling.
% 
% MATLAB R2010b was used for testing
% Need to download the Dimensionality Reduction toolbox (drtoolbox) and
% ensure it's in your path.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


close all;

if ~exist('seed_in', 'var')
    seed_in = [];
end

if ~exist('fileout', 'var') || isempty(fileout)
    fileout = 'output';
end

if ~exist('netw_used', 'var')
    netw_used = [1 3 4 6 7 9 10 11 12 15 16]; % Only draw ROI of interests: numbers correspond to networks
end

if size(netw_used,1) ~= 1
    error('Input argument ''netw_used'' should be a row vector');
end

if ~exist('hemi', 'var')
    hemi = 'lh'; % by default plot left hemisphere only.
end

fileout = char(fileout);

% Dimensionality reduction 
[loc] = Convert2TwoDim(data, seed_in, fileout);

% For visualization purpose, nodes are colored based on their network
% Hence, label of each ROI is needed
[networks, networks_label, lhrh_tag] = Labels_GetNetw;
disp('Drawing graph for :')
disp(char(networks_label(netw_used)))

if(strcmp(hemi, 'rh'))
    lhrh_tag = 3 - lhrh_tag; % DrawGraph only plot lhrh_tag == 1. So set rh ROIs to 1 and lh ROIs to 2.
elseif(strcmp(hemi, 'lhrh'))
    lhrh_tag(:) = 1;
elseif(strcmp(hemi, 'lh'))
    % do nothing
else
    error('hemi not recognized');
end

% Draw graph
lineth = 0.2;
colscale = GrabColorScale;
DrawGraph(data, loc, lineth, networks, colscale, lhrh_tag, netw_used, fileout)



%%%%%%%%%%%%%%%%%%%
% Labels_GetNetw
%%%%%%%%%%%%%%%%%%%
function [Netw, tmplate, lhrh_idx] = Labels_GetNetw

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
    lhrh_tmp = netw{1,2};
    lhrh_tag(i,1) = lhrh_tmp;
    lhrh_netw(i,1) = tmp;
    subnet = netw{1,4};
    if ~isempty(subnet)
        lhrh_subnet(i,1) = subnet;
    end
end


% Get network label according to the template
% Note that Limbic_TempPole and Limbic_OFC are combined as Limbic
tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'Limbic'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'};

for j = 1:numel(tmplate);
    ind = strmatch(tmplate(j), lhrh_netw);
    Netw(ind) = j;
end

hemis = {'LH' 'RH'};
for s = 1:numel(hemis)
    tmp = strmatch(hemis{s}, lhrh_tag);
    lhrh_idx(tmp) = s;
end



%%%%%%%%%%%%%%%%%%%
% Convert2TwoDim
%%%%%%%%%%%%%%%%%%%
function [map] = Convert2TwoDim(data, seed_in, fileout)
% data -    NxN matrix, rows correspond to observations, columns to dimensions
% seed_in - optional, seed for the dimensionality reduction. If not
%           specified, it will be random generation
% fileout - name of output file

if ~exist('seed_in', 'var')
    seed_in = [];
end

% Dimensionality reduction setting
fin_dim = 2;
method = 'SNE';

if (~isempty(seed_in))
    seed_used = seed_in;
else
    seed_used = sum(clock);
end
s = RandStream('mt19937ar','Seed', seed_used);
RandStream.setDefaultStream(s);

% Reduce to 2-dimensional data
[map, ~] = compute_mapping(data, method, fin_dim);

% Save the data (in .mat file) : [x y] coordinates and the seed for
% initialization
mkdir('rawdata')
save(['rawdata/' fileout '.mat'], 'map', 'seed_used');


%%%%%%%%%%%%%%%%%%%
% DrawGraph
%%%%%%%%%%%%%%%%%%%
function DrawGraph(data, loc, lineth, networks, colscale, lhrh_tag, netw_used, fileout)
% data      -  NxN matrix, N refers to number of nodes/observations. Each value
%              representa strength of correlation between nodes. E.g. data(4,5)
%              represents correlation between node 4 and 5
% loc       -  Nx2 matrix, location of nodes on xy plane (from dimensionality reduction)
% lineth    -  threshold of correlation values. Edges with correlation strength
%              below this will not be plotted
% networks  -  Nx1 matrix that represent labels of nodes (in sequence with
%              data). Nodes with same number belong to the same network
% colscale  -  Mx3 matrix, each row represents color for a network, whereby 
%              M = total number of networks (max(networks))
% lhrh_tag  -  Nx1 matrix, optional. Only nodes with value = 1 will be plotted
%              For example, if we only want left hemisphere nodes to be plotted
% netw_used -  1xP, values correspond to network label
%              Example, if only network 4,5,8 need to be plotted, netw_used = [4 5 8];
%              The number refer to 'networks' variable
% fileout   -  name of output file

if size(data,1)~=size(data,2)
    error('Data is not symmetric')
end
if length(loc)~=length(networks)
    error('Different data length found! Unable to map each observation to their network')
end

data = (data+data')./2; % Make sure data is symmetric
data(logical(eye(size(data)))) = 0; 

% Include only the nodes that are of interest, as listed in netw_used  
used_idx = zeros(size(loc,1),1);
for n = netw_used
    used_idx(networks==n) = 1;
end
if ~isempty(lhrh_tag)
    used_idx_lh = [used_idx & lhrh_tag(:)==1]; % only left hemisphere data is used
end

% Draw lines
[data_bin_x data_bin_y] = find(data>lineth);
corr_max = nanmax(data(:));
corr_min = 0;
add_sc = 0; % To adjust the line width, if needed
data = data+add_sc; 
xxx_tot = 0;
xxx = 0;
for s = 1:size(data_bin_x,1)
    xxx_tot = xxx_tot+1;
    if (used_idx_lh(data_bin_x(s))~=0 & used_idx_lh(data_bin_y(s))~=0) % Plot only lines that connect nodes of interest
        if sum(data_bin_y(1:(s-1))==data_bin_x(s) & data_bin_x(1:(s-1))==data_bin_y(s)) == 0 % To make sure only unique lines are plotted
            xxx = xxx+1;
            line([loc(data_bin_x(s),1) loc(data_bin_y(s),1)], [loc(data_bin_x(s),2) loc(data_bin_y(s),2)], 'linewidth', (data(data_bin_x(s), data_bin_y(s))-add_sc), 'color', ([0.3 0.3 0.3] + [0.3 0.3 0.3].*(1-(data(data_bin_x(s), data_bin_y(s))-corr_min)./(corr_max-corr_min))));
            hold on;
        end
    end
end

% Drawing the nodes
for i = netw_used
    used_idx_curr = [networks==i&lhrh_tag==1];
    scatter(loc(used_idx_curr, 1), loc(used_idx_curr, 2), 100, colscale(i,:)./255, 'filled', 'MarkerEdgeColor','k'); hold on
end

% Save the figures
axis equal
axis off
set(gcf, 'PaperPositionMode', 'auto')
mkdir('figures')
print(gcf, '-depsc', '-r600', ['figures/' char(fileout) '.eps']);




%%%%%%%%%%%%%%%%%%%
% GrabColorScale
%%%%%%%%%%%%%%%%%%%
function colorlist = GrabColorScale

colorlist = [205    62    78; ...
   205    62    78; ...
   205    62    78; ...
   205    62    78; ...
   230   148    34; ...
   230   148    34; ...
   230   148    34; ...
   220   248   164; ...
   196    58   250; ...
   196    58   250; ...
     0   118    14; ...
     0   118    14; ...
    70   130   180; ...
    70   130   180; ...
   120    18   134; ...
   120    18   134];

