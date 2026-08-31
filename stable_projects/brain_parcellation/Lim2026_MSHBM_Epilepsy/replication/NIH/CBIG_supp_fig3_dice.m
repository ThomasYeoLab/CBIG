% ====================================
%% Help
% ====================================
%
% CBIG_supp_fig3_dice
%
% This script replicates Supplementary Figure 3: the inter-subject DICE
% coefficient map for the NIH cohort. For each of the 15 networks, pairwise
% DICE coefficients are computed across all NIH subjects that have
% individual-specific parcellations, then averaged. The per-network mean is
% mapped onto the fsaverage6 inflated surface using the group parcellation
% boundaries with a black-red-orange-yellow colormap (range 0.30-0.70).
%
% PREREQUISITES:
%   Individual-specific parcellations must exist in:
%   replication/NIH/input/{subject}/mshbm/NIH/ind_parcellation/test_set/
%     Ind_parcellation_MSHBM_sub1_w80_MRF10.mat
%   Group parcellation: replication/mshbm/NIH/group.mat
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   Individual parcellations: replication/NIH/input/{subject}/mshbm/NIH/
%     ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat
%   Group parcellation: replication/mshbm/NIH/group.mat
%
% OUTPUTS:
%   Prints per-network and mean inter-subject DICE to console.
%   Displays Supplementary Figure 3 on screen (not saved).
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Last tested on 20 May 2026.

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, ...
    'stable_projects/brain_parcellation/Lim2026_MSHBM_Epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');

nih_input_dir = fullfile(data_dir, 'NIH/input');
mshbm_ref_dir = fullfile(proj_dir, 'replication/mshbm/NIH');

addpath(fullfile(CBIG_CODE_DIR, 'utilities/matlab/figure_utilities'));

% ====================================
%% Find subjects with individual parcellations
% ====================================

n_networks     = 15;
parcel_subname = fullfile('mshbm', 'NIH', 'ind_parcellation', 'test_set', ...
    'Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

sub_dirs = dir(fullfile(nih_input_dir, 'sub-*'));
sub_ids  = {};
for s = 1:numel(sub_dirs)
    pfile = fullfile(nih_input_dir, sub_dirs(s).name, parcel_subname);
    if exist(pfile, 'file')
        sub_ids{end+1} = sub_dirs(s).name; %#ok<AGROW>
    end
end
n_subs = numel(sub_ids);
fprintf('\n====== Computing inter-subject DICE (n=%d subjects) ======\n', n_subs);

% ====================================
%% Compute inter-subject DICE
% ====================================

% Load all parcellations
all_labels = cell(n_subs, 1);
for s = 1:n_subs
    pfile = fullfile(nih_input_dir, sub_ids{s}, parcel_subname);
    d     = load(pfile, 'lh_labels', 'rh_labels');
    all_labels{s} = [d.lh_labels(:); d.rh_labels(:)];
end

% Pairwise DICE across all subject pairs
dice_sum = zeros(n_networks, 1);
n_pairs  = 0;

for s1 = 1:n_subs
    for s2 = s1+1:n_subs
        labels1 = all_labels{s1};
        labels2 = all_labels{s2};
        for net = 1:n_networks
            A = (labels1 == net);
            B = (labels2 == net);
            inter = nnz(A & B);
            denom = nnz(A) + nnz(B);
            if denom > 0
                dice_sum(net) = dice_sum(net) + 2 * inter / denom;
            end
        end
        n_pairs = n_pairs + 1;
    end
end

mean_dice_per_net = dice_sum / n_pairs;

fprintf('Inter-subject DICE per network:\n');
for net = 1:n_networks
    fprintf('  Network %2d: %.4f\n', net, mean_dice_per_net(net));
end
fprintf('Mean inter-subject DICE: %.4f\n', mean(mean_dice_per_net));

% ====================================
%% Map DICE onto group parcellation surface
% ====================================

group_data = load(fullfile(mshbm_ref_dir, 'group.mat'), 'lh_labels', 'rh_labels');
lh_labels  = group_data.lh_labels;
rh_labels  = group_data.rh_labels;

lh_DICE = zeros(size(lh_labels));
rh_DICE = zeros(size(rh_labels));
for net = 1:n_networks
    lh_DICE(lh_labels == net) = mean_dice_per_net(net);
    rh_DICE(rh_labels == net) = mean_dice_per_net(net);
end

% ====================================
%% Plot Figure
% ====================================

n    = 256;
n1   = round(n/3);
n2   = round(n/3);
n3   = n - n1 - n2;
r1   = linspace(0, 1, n1)'; g1 = zeros(n1, 1); b1 = zeros(n1, 1);
r2   = ones(n2, 1); g2 = linspace(0, 0.5, n2)'; b2 = zeros(n2, 1);
r3   = ones(n3, 1); g3 = linspace(0.5, 1, n3)'; b3 = zeros(n3, 1);
cmap = [r1 g1 b1; r2 g2 b2; r3 g3 b3];

draw_surface_maps_with_boundary(lh_DICE, rh_DICE, lh_labels, rh_labels, ...
    'fsaverage6', 'inflated', 0.30, 0.70, cmap);
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% ====================================
%% Cleanup
% ====================================

rmpath(fullfile(CBIG_CODE_DIR, 'utilities/matlab/figure_utilities'));

% ====================================
%% Local Functions
% ====================================

function draw_surface_maps_with_boundary(lh_data, rh_data, ...
    lh_labels, rh_labels, mesh_name, surf_type, min_thresh, max_thresh, colors)
% Visualize lh_data and rh_data on the inflated surface with parcel boundaries.
% Boundary vertices are set to min(data) (black). Colorbar uses [min_thresh, max_thresh].

warning('off', 'MATLAB:warn_r14_stucture_assignment');

if ~exist('mesh_name', 'var'); mesh_name = 'fsaverage'; end
if ~exist('surf_type', 'var'); surf_type = 'inflated';  end

pos = [0.020, 0.510, 0.325, 0.470; ...
       0.355, 0.510, 0.325, 0.470; ...
       0.720, 0.760, 0.240, 0.230; ...
       0.720, 0.510, 0.240, 0.230; ...
       0.020, 0.020, 0.325, 0.470; ...
       0.355, 0.020, 0.325, 0.470; ...
       0.720, 0.260, 0.240, 0.230; ...
       0.720, 0.010, 0.240, 0.230];

h = figure; gpos = get(h, 'Position');
gpos(1) = 0; gpos(2) = 0; gpos(3) = 1200; gpos(4) = 600; set(h, 'Position', gpos);

if exist('colors', 'var')
    m = colors / max(colors(:));
    colormap(m);
else
    m = jet;
    colormap(m);
end

for hemis = {'lh' 'rh'}
    hemi = hemis{1};
    mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, surf_type, 'cortex');
    non_cortex = find(mesh.MARS_label == 1);

    if strcmp(hemi, 'lh')
        data   = lh_data;
        labels = lh_labels;
    else
        data   = rh_data;
        labels = rh_labels;
    end

    if size(data,   1) ~= 1; data   = data';   end
    if size(labels, 1) ~= 1; labels = labels'; end

    if size(mesh.vertices, 2) ~= length(data)
        if length(data) == 10242
            from_mesh   = CBIG_ReadNCAvgMesh(hemi, 'fsaverage5', 'sphere', 'cortex');
            target_mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, 'sphere', 'cortex');
            data   = MARS_linearInterpolate(target_mesh.vertices, from_mesh, data);
            labels = MARS_NNInterpolate(target_mesh.vertices, from_mesh, labels);
        else
            error('Not handling %d vertices', length(data));
        end
    end

    if exist('min_thresh', 'var')
        data(data < min_thresh) = min_thresh;
        data(data > max_thresh) = max_thresh;
        data(non_cortex(1)) = min_thresh;
        data(non_cortex(2)) = max_thresh;
    end

    BoundaryVec = zeros(length(labels), 1);
    maxNeighbors = size(mesh.vertexNbors, 1);
    for i = 1:length(labels)
        label_vertex = int32(labels(i));
        for k = 1:maxNeighbors
            v_neighbor = mesh.vertexNbors(k, i);
            if v_neighbor ~= 0 && int32(labels(v_neighbor)) ~= label_vertex
                BoundaryVec(i) = 1;
            end
        end
    end
    data(BoundaryVec == 1) = min(data);

    if strcmp(hemi, 'lh')
        views = {pos(1,:), -90, 0; pos(2,:), 90, 0; pos(3,:), 90, 90; pos(8,:), 90, -90};
    else
        views = {pos(5,:), 90, 0; pos(6,:), -90, 0; pos(4,:), 90, 90; pos(7,:), 90, -90};
    end
    for v = 1:size(views, 1)
        subplot('Position', views{v,1});
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(views{v,2}, views{v,3});
        axis off;
    end
end

if exist('min_thresh', 'var')
    cbax = axes('Position', [0.29 0.5 0.1 0.02], 'visible', 'off');
    caxis(cbax, [min_thresh, max_thresh]);
    colorbar('peer', cbax, 'horiz', 'Position', [0.29 0.5 0.1 0.02]);
end
end

function ncd = revert_shading_interp_behaviour(s, m)
% Revert shading interp behaviour to r2014a style across MATLAB versions.
s    = get(s);
cdat = s.FaceVertexCData;
cl   = get(gca, 'CLim');
sz   = cl(2) - cl(1);
ncd  = zeros(length(cdat), 1, 3);
for x = 1:length(cdat)
    c         = cdat(x, 1);
    idxf      = ((c - cl(1)) / sz) * (size(m, 1) - 1);
    ncd(x,1,:) = m(round(idxf) + 1, :);
end
end
