function [likelihood, results, label, tau] = CBIG_hMRF_update_labels_via_graphcut(full_corr_mat,...
    params, lh_mask, rh_mask, avg_mesh, label, Neighborhood, Smooth_cost, i, tau)
% [likelihood, results, label, tau] = CBIG_hMRF_update_labels_via_graphcut(full_corr_mat,...
% params, lh_mask, rh_mask, avg_mesh, label, Neighborhood, Smooth_cost, i, tau)
%
% This function update parcellation labels via the graphcut algorithm.

% For the notations below:
% N = no of vertices per hemisphere; 
% M = no of cortical vertices for both hemispheres;
% k = total no of parcels 
% 
% Input 
%   - full_corr_mat: (matrix)
%     The 2Nx2N premultiplied matrix computed from multiplying normalized concatenated fMRI data.
%   - params: (struct)
%     The structure containing the input arguments.
%   - lh_mask: (matrix)
%     The Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the left hemisphere.
%   - rh_mask: (matrix)
%     The Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the right hemisphere.
%   - avg_mesh: (struct)
%     The structure containing meshes for both left and right hemispheres.
%   - label: (matrix)
%     The Mx1 resultant parcellation label at the current step.
%   - Neighborhood: (matrix)
%     A MxM sparse matrix containing the neighborhood information for the current mesh.
%   - Smooth_cost: (matrix)
%     A matrix representing smoothcost and to be fed to the GCO algorithm.
%   - i: (integer)
%     An integer controlling the run iterations and termination criteria for GCO.
%     Graphcut will run  10 ^ (i+2) iterations and termination criteria would be tightened to 10^-i% of improvement.
%   - tau: (matrix)
%     An 1xk array representing the current values for the tau hyperparameter.

% Output
%   - likelihood: (matrix)
%     A Mxk matrix representing the current likelihood given a fixed set of labels.
%   - results: (struct)
%     Struct containing some critical intermediate information at the current optimization iteration.
%   - label: (matrix)
%     Resultant Mx1 parcellation label at the current step.
%   - tau: (matrix)
%     An 1xk array representing the current values for the tau hyperparameter. Might differ from initial input tau.

%
% Example
%   - [likelihood, results, label, tau] = CBIG_hMRF_update_labels_via_graphcut(full_corr_mat,...
%      params, lh_mask, rh_mask, avg_mesh, label, Neighborhood, Smooth_cost, 4, tau)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cortex_vertices = sum([lh_mask; rh_mask]);
params.graphCutIterations = 10 ^ (i+2);
termination = 10^-i; % when improvement is less or equal to 10^-i%, the algorithm will terminate the for loop
Enew = 1e18;
convg_log_fileID = fopen(params.convergence_log_path, 'a');

for j = 1:params.num_iterations
    fprintf('Now is iteration %d...\n', j);
    Eold = Enew;
    label = absorb_tiny_clusters(label, params.min_verts_per_cluster, lh_mask, rh_mask, avg_mesh);
    if(params.process_lost_parcels)
        [label, tau] = CBIG_hMRF_find_and_fix_lost_parcels(label, avg_mesh, lh_mask, rh_mask,...
            full_corr_mat, tau, params);
    end
    
    if(params.flip_mismatched_parcels)
        if(rem(j, 2) == 0 || j == params.num_iterations)
            label = detect_flipped_and_change_rh_labels(label, lh_mask, rh_mask,...
                params.num_cluster_per_hemi, params);
        end
    end
    
    [likeli, max_max_likeli, kappa] = compute_global_likelihood_vectorized_kappa(full_corr_mat([lh_mask;...
        rh_mask], [lh_mask; rh_mask]), label', params.num_cluster, params.dim);
    likeli_pos = compute_spatial_likelihood(params, label, lh_mask, rh_mask, tau);
    likeli_sep = compute_label_separation_cost(label, lh_mask, params.num_cluster_per_hemi);
    [label, Enew, E_current_D, E_current_S] = compute_labels_via_graphcut(cortex_vertices, likeli,...
        Neighborhood, Smooth_cost, params, likeli_pos', likeli_sep);
    
    Enew = (Enew - max_max_likeli*cortex_vertices);
    E_current_D = (E_current_D - max_max_likeli*cortex_vertices);
    
    fprintf(convg_log_fileID, 'improvement after %i iterations of %f percent \n', j, (Eold/Enew - 1)*100);
    fprintf(convg_log_fileID, 'Smooth_cost: %f, DataCost: %f, Energy %f \n', E_current_S, E_current_D, Enew);
    
    if ((abs(Eold/Enew-1) * 100) < termination)
        disp('# Hitting termination criteria... #');
        fprintf(convg_log_fileID, 'Hitting termination criteria of %f at iteration %d  \n' , termination, j);
        break
    else
        fprintf(convg_log_fileID, 'Not meeting termination criteria of %f at iteration %d \n' , termination, j);
    end 
end
fclose(convg_log_fileID);

likeli = likeli + max_max_likeli;
likelihood = mean(likeli(label'));
idx = sub2ind(size(likeli), label', 1:max(size(label')));
final_global_likeli = likeli(idx);

results.D = E_current_D;
results.S = E_current_S;
results.kappa = kappa;
results.E = compute_neg_loglikeli(final_global_likeli, results.S);
end

function [likeli, max_max_likeli, kappa] = compute_global_likelihood_vectorized_kappa(PMM_matrix, label, k, dimension) 
% compute the global likelihood term in cost function
miu_times_y = zeros(k, length(label));
for i = 1: k
    indexed_data = find(label == i);
    norm_of_miu = sqrt(sum((sum(PMM_matrix(indexed_data, indexed_data)))));
    miu_times_y(i, :) = sum(PMM_matrix(indexed_data,:), 1);
    miu_times_y(i, :) = miu_times_y(i, :) * 1/norm_of_miu;
    gamma = norm_of_miu / length(indexed_data);
    if(length(indexed_data) == 1)
        if(gamma >= 0.975)
            gamma = 0.975;
        end
    end
    kappa(i) = inv_Ad(dimension, double(gamma));
end

likeli = bsxfun(@plus, compute_log_of_partition_function_vonmises(kappa, dimension)',...
    bsxfun(@times, kappa', miu_times_y));
max_max_likeli = max(max(likeli));
likeli = likeli - max_max_likeli;
end

function [out, exitflag] = inv_Ad(D,gamma)
% Approximating kappa in the global likelihood term
% Please refer to [Lashkari 2010, Discovering structure in the space of fMRI selectivity profiles]
% for the detailed steps for approximating kappa

outu = (D-1)*gamma/(1-gamma^2) + D/(D-1)*gamma; % note that we already reduced D by 1
normalization_z = besseli(D/2-1, outu);

if ((normalization_z == Inf)||(isnan(normalization_z))||(normalization_z == 0)) % 
    out = outu - D/(D-1)*gamma/2;
    exitflag = Inf;
else
    options = optimset('Display','off');
    [outNew, ~, exitflag] = fzero(@(argum) besseli_Ad(argum,D)-gamma, outu, options);
    if exitflag == 1
        out = outNew;
    else
        out = outu - D/(D-1)*gamma/2;
    end
end
end

function out = besseli_Ad(in,D)
% The Ad function, refer to [Lashkari 2010, Discovering structure in the space of fMRI selectivity profiles]
out = besseli(D/2,in) ./ besseli(D/2-1,in);
end

function likeli = compute_spatial_likelihood(params, label, lh_mask, rh_mask, tau)
% compute the spatial likelihood term in cost function

[lh_sphere_xyz, rh_sphere_xyz] = get_spherical_coordinates(params.mesh_type);
lh_sphere_xyz = lh_sphere_xyz(:, lh_mask == 1)';
rh_sphere_xyz = rh_sphere_xyz(:, rh_mask == 1)';
lh_sphere_xyz = bsxfun(@rdivide, lh_sphere_xyz, sqrt(sum(lh_sphere_xyz.^2, 2)));
rh_sphere_xyz = bsxfun(@rdivide, rh_sphere_xyz, sqrt(sum(rh_sphere_xyz.^2, 2)));
xyz_dist = [lh_sphere_xyz; rh_sphere_xyz];

if(size(label,1) > 1)
    label = label';
end

n = length(label);
t1 = 1: 1: n;
t2 = label;

vertex_cluster_match = zeros(n, params.num_cluster);
vertex_cluster_match(sub2ind(size(vertex_cluster_match), t1, t2)) = 1;

mean_spatial_dir = bsxfun(@times, (xyz_dist'*vertex_cluster_match),...
    1./sqrt(sum((xyz_dist'*vertex_cluster_match).^2)));
log_norm_partition_func(tau==0) = 0;
log_norm_partition_func(tau~=0) = compute_log_of_partition_function_vonmises(single(tau(tau~=0)), 3);

% compute likelihood c(tau) + tau*nu'*x'
likeli_xyz_dis = bsxfun(@plus, log_norm_partition_func, bsxfun(@times, tau, (mean_spatial_dir'*xyz_dist')'));

max_max_likeli_xyz = max(max(likeli_xyz_dis));
likeli_xyz_dis = likeli_xyz_dis - max_max_likeli_xyz;
likeli = params.w_xyz * likeli_xyz_dis;
end

function [lh_sphere_xyz, rh_sphere_xyz] = get_spherical_coordinates(mesh_type)
% get the coordinates on the sphere of a given surface mesh
if(contains(mesh_type, 'fs_LR'))
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh_type, 'sphere');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh_type, 'sphere');
else
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_type, 'sphere', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_type, 'sphere', 'cortex');
end

lh_sphere_xyz = lh_avg_mesh.vertices;
rh_sphere_xyz = rh_avg_mesh.vertices;
end


function cost = compute_label_separation_cost(label, lh_mask, num_cluster_per_hemi)
% compute the label separation cost term
lh_cortex_vertices = sum(lh_mask);
cost = zeros(num_cluster_per_hemi*2, length(label));

cost(1:num_cluster_per_hemi, 1:lh_cortex_vertices) = 1e10;
cost(num_cluster_per_hemi+1:end, 1:lh_cortex_vertices) = 0;

cost(1:num_cluster_per_hemi, lh_cortex_vertices+1:end) = 0;
cost(num_cluster_per_hemi+1:end, lh_cortex_vertices+1:end) = 1e10;

max_likeli = max(cost(:));
cost = cost - max_likeli;
end

function [label, Enew, D, S] = compute_labels_via_graphcut(cortex_vertices, probs, Neighborhood,...
    Smooth_cost, params, likeli_pos, likeli_sep)
% compute labels based on graphcut algorithm

h = GCO_Create(cortex_vertices, params.num_cluster);
normalize = cortex_vertices;

data_cost = single(-(1/normalize) * ((probs+eps))+eps);
pos_cost = single(-(1/normalize) * ((likeli_pos+eps))+eps);
sep_cost = single(-(1/normalize) * (likeli_sep + eps) + eps);
GCO_SetDataCost(h, data_cost+pos_cost+sep_cost);
GCO_SetNeighbors(h, Neighborhood);
GCO_SetSmoothCost(h, Smooth_cost * single((1/normalize)));
GCO_Expansion(h, params.graphCutIterations);
label = GCO_GetLabeling(h);
[Enew, D, S] = GCO_ComputeEnergy(h);

GCO_Delete(h);
Enew = Enew * normalize;
D = D * normalize;
S = S * normalize;
end

function negative_log_likeli = compute_neg_loglikeli(global_likeli, smoothcost)
% compute negative log likelihood
negative_log_likeli = - sum(global_likeli) + smoothcost;
end

function label = detect_flipped_and_change_rh_labels(label, lh_mask, rh_mask, num_cluster_per_hemi, params)
% flipped parcels refer to those with swapped spatical locations that disrupt the desired homotopy
% e.g., L1 and R1, L2 and R2 are supposed to be homotopic pairs but instead, L1 is matched to R2 and L2 to R1
% we use this function to detect and fix these parcels.

disp('##################################');
disp('# Flipping mismatched parcels... #');
disp('##################################');

[lh_full_label, rh_full_label] = CBIG_hMRF_get_left_right_overlapping_labels(lh_mask, rh_mask,...
    label, num_cluster_per_hemi);

% detect flipped pairs
assign = hungarian_cluster_match(lh_full_label, rh_full_label, params);
pair_to_flip = find_parcel_to_flip(assign, rh_full_label, rh_full_label);
rh_full_flipped = flip_rh_parcels_by_lh(rh_full_label, pair_to_flip);

% restore rh label such that it does not overlap with lh label
rh_full_flipped(rh_full_flipped ~= 0) = rh_full_flipped(rh_full_flipped ~= 0) + num_cluster_per_hemi;

% update label vector
label = [lh_full_label(lh_mask) rh_full_flipped(rh_mask)];
end

function pair_to_flip = find_parcel_to_flip(assign, lh_full_label, rh_full_label)
% use hungarian matching to flip parcel pairs

pair_to_flip = [];
if_detected = zeros(1,max(lh_full_label)); % to ensure that we don't pick the circular pairs multiple times
for lh_cluster_idx = 1:max(lh_full_label)
    rh_matched_cluster_idx = assign(lh_cluster_idx);
    if(rh_matched_cluster_idx ~= lh_cluster_idx)
        if(assign(rh_matched_cluster_idx) == lh_cluster_idx && ~if_detected(lh_cluster_idx)) % circular mismatch
            % check if "change in total correlation does not drop by ?%" criteria
            cluster_A_idx = lh_cluster_idx;
            cluster_B_idx = rh_matched_cluster_idx;
            lh_A_verts = lh_full_label == cluster_A_idx;
            rh_A_verts = rh_full_label == cluster_A_idx;
            lh_B_verts = lh_full_label == cluster_B_idx;
            rh_B_verts = rh_full_label == cluster_B_idx;
            
            pair_to_flip = [pair_to_flip; lh_cluster_idx, assign(lh_cluster_idx)];

            if_detected(lh_cluster_idx) = 1;
            if_detected(assign(lh_cluster_idx)) = 1;
        end
    end
end
end

function assign = hungarian_cluster_match(lh_full_label, rh_full_label, params)
% use Hungarian/Munkres algorithm to assign each rh parcel to lh based on maximum vertex overlap
load(params.lh_vert2rh_vert, 'lh_one2rh_one');
cost_matrix = zeros(max(lh_full_label), max(rh_full_label));
for row_cluster_idx = 1:max(lh_full_label)
    for col_cluster_idx = 1:max(rh_full_label)
        % count number of overlapping vertices for the current cluster
        row_cluster_mapped_to_rh_verts = lh_one2rh_one(lh_full_label == row_cluster_idx);
        col_cluster_rh_verts = find(rh_full_label == col_cluster_idx);
        cost_matrix(row_cluster_idx, col_cluster_idx) = -length(intersect(row_cluster_mapped_to_rh_verts,...
            col_cluster_rh_verts));
    end
end
[assign, ~] = munkres(cost_matrix);
end

function rh_full_label_flipped = flip_rh_parcels_by_lh(rh_full_label, pair_to_flip)
% flip the parcels on right hemispheres to match left hemisphere topology
rh_full_label_flipped = rh_full_label;
for cur_pair_idx = 1:size(pair_to_flip,1)
    cluster_A = pair_to_flip(cur_pair_idx,1);
    cluster_B = pair_to_flip(cur_pair_idx,2);
    %fix lh, flip rh
    rh_full_label_flipped(rh_full_label== cluster_A) = cluster_B;
    rh_full_label_flipped(rh_full_label== cluster_B) = cluster_A;
end
end

function label = absorb_tiny_clusters(label, min_verts_per_cluster, lh_mask, rh_mask, avg_mesh)
% looks for components of <= min_verts_per_cluster and merge these components
disp('###########################');
disp('# Absorb tiny clusters... #');
disp('###########################');

lh_avg_mesh = avg_mesh.lh_avg_mesh;
rh_avg_mesh = avg_mesh.rh_avg_mesh;
lh_label = label(1: sum(lh_mask));
rh_label = label(sum(lh_mask)+1: end);

lh_full_label(lh_mask) = lh_label;
rh_full_label(rh_mask) = rh_label;

lh_Nbors = lh_avg_mesh.vertexNbors;
rh_Nbors = rh_avg_mesh.vertexNbors;

lh_full_label = reassign_small_cluster_label(lh_full_label, lh_avg_mesh, min_verts_per_cluster, lh_Nbors);
rh_full_label = reassign_small_cluster_label(rh_full_label, rh_avg_mesh, min_verts_per_cluster, rh_Nbors);

label = [lh_full_label(lh_mask)';rh_full_label(rh_mask)'];

end

function hemi_full_label = reassign_small_cluster_label(hemi_full_label, hemi_avg_mesh, min_verts_per_cluster,...
    hemi_nbors)
% This function picks the components that consist of vertices less or equal to min_verts_per_cluster, and then
% reassign the component cluster labels by referring to the neighbors of the current clique.
exist_reassigned_clique = 1;
while(exist_reassigned_clique) 
    [hemi_ci, ~] = CBIG_hMRF_generate_components_one_hemi(hemi_avg_mesh, hemi_full_label);
    [exist_reassigned_clique, hemi_full_label] = search_for_small_clique_until_one_is_reassigned(hemi_ci,...
        min_verts_per_cluster, hemi_full_label, hemi_nbors);
end
end

function [exist_reassigned_clique, hemi_full_label] = search_for_small_clique_until_one_is_reassigned(hemi_ci,...
min_verts_per_cluster, hemi_full_label, hemi_nbors)
% goes through all connected components and try to merge small ones, return when a merging is complete

exist_reassigned_clique = 0;

for i = 1: length(unique(hemi_ci))
    vertices_in_cur_clique = find(hemi_ci == i);
    if(length(vertices_in_cur_clique) <= min_verts_per_cluster)
        [hemi_full_label, cur_clique_reassigned] = update_label_for_one_clique(hemi_full_label,...
            vertices_in_cur_clique, hemi_nbors);
        if(cur_clique_reassigned)
            exist_reassigned_clique = 1;
            break;
        end
    end
end
end

function [hemi_full_label, clique_reassigned] = update_label_for_one_clique(hemi_full_label,...
vertices_in_cur_clique, hemi_nbors)
% reassign the "clique" labels by referring to the neighbors of the current clique
all_neighbors_of_clique = [];
for j = 1: length(vertices_in_cur_clique)
    all_neighbors_of_vertex = nonzeros(hemi_nbors(:, vertices_in_cur_clique(j)));
    cluster_labels_of_neighbors = hemi_full_label(all_neighbors_of_vertex);
    all_neighbors_of_vertex(cluster_labels_of_neighbors == 0) = []; 
    all_neighbors_of_clique = [all_neighbors_of_clique; all_neighbors_of_vertex];
end
all_valid_neighbors_of_clique = unique(all_neighbors_of_clique);
all_valid_neighbors_of_clique = setdiff(all_valid_neighbors_of_clique, vertices_in_cur_clique);
cluster_labels_of_all_valid_neighbors = hemi_full_label(all_valid_neighbors_of_clique);

if(~isempty(cluster_labels_of_all_valid_neighbors))
    hemi_full_label(vertices_in_cur_clique) = mode(cluster_labels_of_all_valid_neighbors);
    clique_reassigned = 1;
else
    clique_reassigned = 0;
end

end

function log_out = compute_log_of_partition_function_vonmises(k, d)
% Computes the logarithm of the partition function of vonMises-Fisher as a function of k
% k is the concentration parameter of the vonmises function. d is the dimensionality.

sizek = size(k);
k = k(:);
log_out = (d/2-1).*log(k) - log(besseli((d/2-1)*ones(size(k)),k));
if (d <= 3) % for the case of xyz term
    k0 = 25;
else
    k0 = CBIG_hMRF_initialize_lambda_in_vonmises_partition_func(0,d);
end

fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
if isinf(fk0)
    k0 = 0.3313 * d;
    fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
end
nGrids = 1000;

maskof = find(max((k>k0),(isinf(log_out))));
nkof = length(maskof);

if (nkof > 0)
    kof = k(maskof);
    ofintv = (kof - k0)/nGrids;
    tempcnt = (1:nGrids) - 0.5;
    ks = k0 + repmat(tempcnt,nkof,1).*repmat(ofintv,1,nGrids);
    CBIG_Adsum = sum(1./((0.5*(d-1)./ks) + sqrt(1+(0.5*(d-1)./ks).^2)) ,2);
    log_out(maskof) =  fk0 - ofintv .* CBIG_Adsum;
end

log_out = reshape(log_out,sizek);
end
