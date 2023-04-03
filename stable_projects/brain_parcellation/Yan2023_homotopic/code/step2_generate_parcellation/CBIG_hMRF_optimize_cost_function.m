function results = CBIG_hMRF_optimize_cost_function(full_corr_mat, dim, params)
% results = CBIG_hMRF_optimize_cost_function(full_corr_mat, dim, params)
%
% This is the main clustering algorithm of generating the parcellation.

% For the notations below:
% N = no of vertices per hemisphere; 

% Input
%   - full_corr_mat (matrix):
%     The 2Nx2N input premultiplied fMRI matrix.
%   - dim (double):
%     The total number of time points of concatenated fMRI time courses, before taking the inner product.
%   - params: (struct)
%     The structure containing the input arguments.
% Output
%   - results: (struct)
%     A struct containing the final parcellation labels and various additional results
%
% Example
%   - results = CBIG_hMRF_optimize_cost_function(full_corr_mat, 300000, params)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

params.dim = dim - 1; % total length of concatenated fmri data; -1 to account for lashkari(2010) approximation

disp('###############################');
disp('# Parameter initialization... #');
disp('###############################');
[lh_avg_mesh, rh_avg_mesh, lh_mask, rh_mask] = CBIG_hMRF_load_mesh_mask_by_mesh_type(params.mesh_type);
avg_mesh.lh_avg_mesh = lh_avg_mesh;
avg_mesh.rh_avg_mesh = rh_avg_mesh;
avg_mask = [lh_mask; rh_mask];

disp('Form the neighborhood matrix...');
[Neighborhood_lh, Neighborhood_rh, Neighborhood_lh_rh] = initialize_neighborhood_components(params,...
lh_avg_mesh, rh_avg_mesh, lh_mask, rh_mask);
Neighborhood = CBIG_hMRF_update_whole_brain_neighborhood(Neighborhood_lh, Neighborhood_rh, Neighborhood_lh_rh,...
    params.c, params.d);

disp('Form the smoothcost matrix...');
smooth_cost_mat = CBIG_hMRF_initialize_homotopic_smoothcost_mat(params.num_cluster_per_hemi, params.d);

disp('Load the border matrix as cluster initialization input...');
grad_mat = load_gradient_matrix_from_file(params, lh_mask, rh_mask);

% fixate the random initialization for each params.seed for initialization stage
rng(params.seed, 'twister'); 
disp('Initializing tau as a uniform vector...');
tau = zeros(1,params.num_cluster) + CBIG_hMRF_initialize_lambda_in_vonmises_partition_func(params.initial_tau,3);

disp('##########################');
disp('# Cluster initialization #');
disp('##########################');

disp('Random initialization over the cortex using lh gradient data...');
[~, rh_likeli] = initialize_global_likelihood(params.num_cluster_per_hemi, params.dim,...
    full_corr_mat(avg_mask, avg_mask), grad_mat, sum(lh_mask));
disp('Assigning labels symmetrically over the 2 hemispheres...');
label = initialize_label_on_rh_then_mirror_to_lh(rh_likeli, Neighborhood, params, avg_mask,...
    avg_mesh, full_corr_mat);

%% ensure initialization contains full label on either hemi
[label, tau] = CBIG_hMRF_find_and_fix_lost_parcels(label, avg_mesh, lh_mask, rh_mask, full_corr_mat, tau, params);
initial_label = label;

for i = [0 2 4]
    [likelihood, results, label, tau] = CBIG_hMRF_update_labels_via_graphcut(full_corr_mat, params,...
    lh_mask, rh_mask, avg_mesh, label, Neighborhood, smooth_cost_mat, i, tau);
end

results.initial_full_label = zeros(1, max(size(avg_mask)));
results.initial_full_label(avg_mask == 1) = initial_label(1: (sum(lh_mask) + sum(rh_mask)));
save_intermediate_results(label, params, lh_mask, rh_mask, results, 'pre_hyperpara_autotune');

if(params.decrease_d)
    disp('########################################');
    disp('# relax LR symmetry by decreasing d... #');
    disp('########################################');
    initial_d = params.d;
    [likelihood, results, label, params, tau] = CBIG_hMRF_adapt_d(params, label, 4, full_corr_mat,...
        lh_mask, rh_mask, avg_mesh, tau, Neighborhood_lh, Neighborhood_rh, Neighborhood_lh_rh);
    final_d = params.d;
    save_intermediate_results(label, params, lh_mask, rh_mask, results, 'after_decrease_d');
end

if(params.decrease_c)
    disp('##############################################');
    disp('# relax grad-weighted MRF by decreasing c... #');
    disp('##############################################');
    initial_c = params.c;
    [likelihood, results, label, params, tau] = CBIG_hMRF_adapt_c(params, label, 4, full_corr_mat,...
        lh_mask, rh_mask, avg_mesh, smooth_cost_mat, tau, Neighborhood_lh, Neighborhood_rh, Neighborhood_lh_rh);
    save_intermediate_results(label, params, lh_mask, rh_mask, results, 'after_decrease_c');
    final_c = params.c;
end

if(params.increase_tau)
    disp('##################################################');
    disp('# Perform increase tau for parcels that split... #');
    disp('##################################################');
    [~, results, label] = CBIG_hMRF_adapt_tau(params, label, full_corr_mat, lh_mask, rh_mask,...
    avg_mesh, Neighborhood, smooth_cost_mat, tau, 4, likelihood, results);
    save_intermediate_results(label, params, lh_mask, rh_mask, results, 'after_increase_tau');
end

%% write to results struct
results.full_label(avg_mask) = label;
results.lh_label = results.full_label(1:length(lh_mask));
results.rh_label = results.full_label(length(lh_mask)+1:end);

if(params.decrease_d)
    results.initial_d = initial_d;
    results.final_d = final_d;
end

if(params.decrease_c)
    results.initial_c = initial_c;
    results.final_c = final_c;
end

fprintf('Result for seed %d has been generated.\n', params.seed);
end

function label = initialize_label_on_rh_then_mirror_to_lh(probs, Neighborhood, params,...
    avg_mask, avg_mesh, full_corr_mat)
% initialize parcellation label on the right hemisphere then reflect to the left

%% initialize on rh
lh_mask = avg_mask(1: length(avg_mask)/2);
rh_mask = avg_mask(length(avg_mask)/2 + 1: end);

num_cortical_verts = length(rh_mask);
lh_cortex_vertices = sum(sum(lh_mask));
rh_cortex_vertices = sum(sum(rh_mask));
Neighborhood_rh = Neighborhood(lh_cortex_vertices+1:end, lh_cortex_vertices+1:end);

h = GCO_Create(rh_cortex_vertices, params.num_cluster_per_hemi);
data_cost = single(-(1/rh_cortex_vertices) * ((probs(:,lh_cortex_vertices+1:end)+eps))+eps);

GCO_SetDataCost(h, data_cost);

Smooth_cost = ones(params.num_cluster_per_hemi, params.num_cluster_per_hemi, 'single');
Smooth_cost(1: (params.num_cluster_per_hemi+1): end) = 0;

GCO_SetNeighbors(h, Neighborhood_rh);
GCO_SetSmoothCost(h, Smooth_cost * single(params.d * (1/rh_cortex_vertices)));
GCO_Expansion(h, params.graphCutIterations);
rh_label = GCO_GetLabeling(h);
GCO_Delete(h);

%% mirror rh initialization to lh
lh_avg_mesh = avg_mesh.lh_avg_mesh;
rh_avg_mesh = avg_mesh.rh_avg_mesh;
if(contains(params.mesh_type, 'fs_LR'))
    [lh_full_label, rh_full_label] = mirror_rh_label_to_lh_for_initialization_fsLR(lh_mask, rh_mask,...
        rh_label, lh_avg_mesh);
else
    [lh_full_label, rh_full_label] = mirror_rh_label_to_lh_for_initialization_fs6(lh_mask, rh_mask,...
        rh_label, lh_avg_mesh, params);
end

%% merge distant parcel pairs
distant_pairs = find_spatially_distant_pairs(lh_full_label, rh_full_label, params.num_cluster_per_hemi,...
    params.mesh_type, 50, params.output_folder);
lh_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(full_corr_mat(1:num_cortical_verts,...
    1:num_cortical_verts), distant_pairs, lh_full_label, lh_avg_mesh, lh_mask);
rh_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(full_corr_mat(num_cortical_verts+1:end,...
    num_cortical_verts+1:end), distant_pairs, rh_full_label, rh_avg_mesh, rh_mask);
lh_label = lh_full_label(lh_mask);

rh_label = int32(rh_full_label(rh_mask)) + params.num_cluster_per_hemi;
label = [lh_label rh_label];
end

function [lh_full_label, rh_full_label] = mirror_rh_label_to_lh_for_initialization_fsLR(lh_mask, rh_mask,...
    rh_label, lh_avg_mesh)
% this function is for mirroring initial rh labels to lh labels for fs_LR_32k surface.

rh_full_label(rh_mask) = rh_label;

medialwall_union = (rh_mask & lh_mask);
lh_full_label(medialwall_union) = rh_full_label(medialwall_union);

rh_not_lh_medial_verts = (rh_mask + (rh_mask | lh_mask));
% 0 - medial wall on both hemi; 1 - medial wall on rh only; 2 - not medial on either hemi
rh_not_lh_medial_verts = find(rh_not_lh_medial_verts == 1);
remaining_rh_not_lh_medial_verts = rh_not_lh_medial_verts;

% some patches of vertices all share a label of zero; the assignment would start from the boundary then move inward
while(~isempty(remaining_rh_not_lh_medial_verts))

    medial_verts_attempt_to_assign = remaining_rh_not_lh_medial_verts;

    for i = 1: length(medial_verts_attempt_to_assign)
        vert_idx = medial_verts_attempt_to_assign(i);
        idx_Nbors = nonzeros(lh_avg_mesh.vertexNbors(:, vert_idx));
        val_Nbors = nonzeros(lh_full_label(idx_Nbors));
        if(~isempty(val_Nbors))
            lh_full_label(vert_idx) = mode(val_Nbors);
            remaining_rh_not_lh_medial_verts(remaining_rh_not_lh_medial_verts == vert_idx) = [];
        end
    end
end
end

function [lh_full_label, rh_full_label] = mirror_rh_label_to_lh_for_initialization_fs6(lh_mask, rh_mask,...
    rh_label, lh_avg_mesh, params)
% this function is for mirroring initial rh labels to lh labels for fsaverage6 surface.

fake_label_to_remove = 10000;
lh_full_label = ones(length(lh_mask),1) * fake_label_to_remove;
rh_full_label(rh_mask) = rh_label;

% map from rh labels to lh labels
if(exist('params','var'))
    load(params.rh_vert2lh_vert, 'rh_one2lh_one');
end

for rh_vert_idx = 1:length(rh_mask)
    if(rh_one2lh_one(rh_vert_idx) && rh_full_label(rh_vert_idx) ~= 0) 
        % do the assignment only if cur rh vert is successfully matched to a lh vert
        % AND the current rh vert should NOT be a medial wall vert
        lh_full_label(rh_one2lh_one(rh_vert_idx)) = rh_full_label(rh_vert_idx);
    end
end
lh_full_label = assign_empty_verts_initialization(lh_full_label, lh_avg_mesh, fake_label_to_remove);
lh_full_label(~lh_mask) = 0;
end

function hemi_label_merged = assign_empty_verts_initialization(hemi_full_label, hemi_avg_mesh, fake_label_to_remove) 
% assign empty vertices in initialized label on a given hemisphere

if (max(size(hemi_full_label)) ~= max(size(hemi_avg_mesh.vertices)))
    error('number of vertices in label vector dont match mesh vertices');
end

if(size(hemi_full_label, 1) ~= 1)
    warning('assume hemi_full_label is 1 x N rotating');
    hemi_full_label = hemi_full_label';
end

hemi_label_merged = hemi_full_label;
hemi_label_to_merge_indices = (hemi_label_merged == fake_label_to_remove);

while(max(hemi_label_to_merge_indices))
    % need to run merge_vertex_label for multiple rounds until undesired labels are completely merged
    hemi_label_merged = merge_vertex_label(hemi_label_merged, hemi_label_to_merge_indices, hemi_avg_mesh);
    hemi_label_to_merge_indices = (hemi_label_merged == fake_label_to_remove);
end
end

function hemi_label_merged = merge_vertex_label(hemi_full_label, hemi_label_to_merge_indices, hemi_avg_mesh)
% `hemi_label_to_merge_indices` indicates the positions in `hemi_full_label` where the labels need to be merged
% by `merge` it means to remove these vertices one by one, i.e., to replace the label of such a vertex
% with the a most frequent and different label of its neighbors

hemi_idx = find(hemi_label_to_merge_indices > 0);
hemi_comp_idx = find(hemi_label_to_merge_indices == 0); 
hemi_label_merged = hemi_full_label;

for i = 1:length(hemi_idx)
    current = hemi_idx(i);
    nbors = hemi_avg_mesh.vertexNbors(:,current);
    nbors_intersect = intersect(nbors,hemi_comp_idx);
    value_to_assign = mode(hemi_full_label(nbors_intersect));
    if(~isempty(nbors_intersect)) % if this vertex is currently within a region that needs to be reassigend
        if(isempty(intersect(value_to_assign,hemi_full_label(nbors_intersect))))
            value_to_assign = hemi_full_label(nbors_intersect(1));
        end
    else
        value_to_assign = hemi_full_label(current);
    end
    hemi_label_merged(current) = value_to_assign;
end
end

function distant_pairs = find_spatially_distant_pairs(lh_full_label, rh_full_label,...
    num_cluster_per_hemi, mesh_type, dist_threshold, output_folder)
% find the spatially distant homotopic pairs for a given parcellation.
% since the fsaverage6 surface mesh is by itself asymmetric, it is easier to compute and compare the coordinates 
% parcel centre by first projecting the given fsaverage6 parcellation labels to the fs_LR_32k surface.

if(size(lh_full_label, 1) ~= 1)
    lh_full_label = lh_full_label';
end

if(size(rh_full_label, 1) ~= 1)
    rh_full_label = rh_full_label'; 
end

if(max(rh_full_label) > max(lh_full_label))
    rh_full_label = rh_full_label - num_cluster_per_hemi;
end

if(strcmp(mesh_type, 'fs_LR_32k'))
    ;
elseif(strcmp(mesh_type, 'fsaverage6'))
    folder_to_write = fullfile(output_folder, ['temp_folder' num2str(rand(1))]);
    [lh_full_label, rh_full_label] = CBIG_project_fsaverage2fsLR(lh_full_label, rh_full_label,...
        mesh_type, 'label', folder_to_write);
else
    error('Input mesh type not recognized.');
end


num_pair = 0;
metric = ones(num_cluster_per_hemi, 1) * -1;
% by default it's -1; unless the current pair exist, we replace -1 with the within pair distance
for cur_parcel_idx = 1: num_cluster_per_hemi
    if(sum(lh_full_label == cur_parcel_idx) ~= 0)
        if(sum(rh_full_label == cur_parcel_idx) ~= 0)
            num_pair = num_pair + 1;
            [lh_index] = find(lh_full_label == cur_parcel_idx);
            [rh_index] = find(rh_full_label == cur_parcel_idx);
            overlap = intersect(lh_index, rh_index);
            metric(cur_parcel_idx) = (2 * length(overlap) / length([lh_index; rh_index]));
        end
    end
end

lh_avg_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k', 'sphere');
rh_avg_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k', 'sphere');
lh_spatial_cor = lh_avg_mesh.vertices;
rh_spatial_cor = rh_avg_mesh.vertices;

distant_pairs = [];
for cur_parcel_idx = 1:num_cluster_per_hemi
    if(metric(cur_parcel_idx) == 0)
        fprintf('There is no overlapping vertices for the current parcel pair.\n');
        lh_cur_parcel_centre = abs(mean(lh_spatial_cor(:, lh_full_label == cur_parcel_idx),2));
        rh_cur_parcel_centre = abs(mean(rh_spatial_cor(:, rh_full_label == cur_parcel_idx),2));
        distance = norm(lh_cur_parcel_centre - rh_cur_parcel_centre);
        
        if(distance > dist_threshold)
            distant_pairs = [distant_pairs; cur_parcel_idx];
        end
    end
end
end

function [lh_likeli, rh_likeli] = initialize_global_likelihood(num_cluster_per_hemi, dimension, PMM_matrix,...
    grad_matrix, lh_cortex_vertices) 
% to initialize global likelihood on both hemispheres based on given PMM and gradient matrices

% using randomly selected vertices for initialization
grad_matrix = mean(grad_matrix, 1);
lh_grad_matrix = grad_matrix(1: lh_cortex_vertices);
rh_grad_matrix = grad_matrix(lh_cortex_vertices+1:end);
lh_low_grad_idx = find(lh_grad_matrix < 0.05);
rh_low_grad_idx = find(rh_grad_matrix < 0.05) + lh_cortex_vertices;
kappa = CBIG_hMRF_initialize_lambda_in_vonmises_partition_func(0, dimension);

% choose num_cluster_per_hemi verts out of the low gradient verts
lh_rand_k_verts = lh_low_grad_idx(datasample(1: length(lh_low_grad_idx), num_cluster_per_hemi, 'Replace', false));
lh_miu_times_x = PMM_matrix(lh_rand_k_verts,:);
lh_likeli = single(kappa * lh_miu_times_x);
lh_max_max_likeli = max(max(lh_likeli));
lh_likeli = lh_likeli - lh_max_max_likeli;

% choose num_cluster_per_hemi verts out of the low gradient verts
rh_rand_k_verts = rh_low_grad_idx(datasample(1: length(rh_low_grad_idx), num_cluster_per_hemi, 'Replace', false));
rh_miu_times_x = PMM_matrix(rh_rand_k_verts,:);
rh_likeli = single(kappa * rh_miu_times_x);
rh_max_max_likeli = max(max(rh_likeli));
rh_likeli = rh_likeli - rh_max_max_likeli;
end

function grad_mat = load_gradient_matrix_from_file(params, lh_mask, rh_mask)
% to load gradient matrices from given file directory

lh_grad_struct = load(params.lh_grad_file);
rh_grad_struct = load(params.rh_grad_file);
lh_grad_matrix = lh_grad_struct.border_matrix(:, lh_mask);
rh_grad_matrix = rh_grad_struct.border_matrix(:, rh_mask);
grad_mat = [lh_grad_matrix, rh_grad_matrix];
end

function [Neighborhood_lh, Neighborhood_rh, Neighborhood_lh_rh] = initialize_neighborhood_components(params, ...
lh_avg_mesh, rh_avg_mesh, lh_mask, rh_mask)
% initialize 3 neighborhood components: Neighborhood_lh, Neighborhood_rh, Neighborhood_lh_rh

% the setting used in the paper with gordon gradient prior fused in the MRF
load(params.lh_grad_file, 'border_matrix');
Neighborhood_lh = build_sparse_gradient(lh_avg_mesh, lh_mask, border_matrix, params.k);
load(params.rh_grad_file, 'border_matrix');
Neighborhood_rh = build_sparse_gradient(rh_avg_mesh, rh_mask, border_matrix, params.k);

if(strcmp(params.mesh_type, 'fs_LR_32k'))
    Neighborhood_lh_rh = diag(ones(length(lh_avg_mesh.vertexNbors),1)); % connecting homotopic vertices
    Neighborhood_lh_rh(~lh_mask, :) = [];
    Neighborhood_lh_rh(:, ~rh_mask) = [];
else
    load(params.nborhood_file, 'Neighborhood_lh_rh');
    Neighborhood_lh_rh(~lh_mask, :) = [];
    Neighborhood_lh_rh(:, ~rh_mask) = [];
end
end

function sparse_gradient = build_sparse_gradient(avg_mesh, avg_mask, border_matrix, k) 
% to build sparse representation of input gradient matrix, taking the exponential k into account (see supplementary)

vertices = max(size(avg_mesh.vertexNbors));
r = reshape(repmat(13: vertices, 6, 1), 1, 6*(vertices-12));
sparse_gradient = sparse(double(r), double(reshape(avg_mesh.vertexNbors(1:6, 13:end),...
[1, 6*(vertices-12)])), reshape(compute_stable_E(border_matrix(:,13: end), k), (size(r))));

for i = 1:12
    sparse_gradient(i, avg_mesh.vertexNbors(1:5, i)) = compute_stable_E(border_matrix(1:5, i), k);
end

cortical_vertices = find(avg_mask == 1);
sparse_gradient = sparse_gradient(cortical_vertices, cortical_vertices);
end

function x = compute_stable_E(x, k)
% computes a stable version of given x, since >1 or <0 can easily incur numeric overflow in the exponential step
x(x > 1) = 1;
x(x < 0) = 0;
x = exp(-k*x) - exp(-k);
end

function save_intermediate_results(label, params, lh_mask, rh_mask, results, intermediate_step_name)
% to save intermediate parcellation results

lh_full_label(lh_mask) = label(1:sum(lh_mask));
rh_full_label(rh_mask) = label(sum(lh_mask)+1:end);

results.full_label = [lh_full_label rh_full_label];
results.lh_label = lh_full_label;
results.rh_label = rh_full_label;

save(fullfile(params.debug_out_folder, [params.output_name '_seed_' num2str(params.seed) '_'...
    intermediate_step_name '.mat']), 'results', 'params');
end