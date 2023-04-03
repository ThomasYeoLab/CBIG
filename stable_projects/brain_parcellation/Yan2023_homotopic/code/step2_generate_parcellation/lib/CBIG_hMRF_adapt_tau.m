function [likelihood, results, label] = CBIG_hMRF_adapt_tau(params, label, full_corr_mat, lh_mask, rh_mask,...
    avg_mesh, Neighborhood, Smooth_cost, tau, i, likelihood, results)
% [likelihood, results, label] = CBIG_hMRF_adapt_tau(params, label, full_corr_mat, lh_mask, rh_mask,...
% avg_mesh, Neighborhood, Smooth_cost, tau, i, likelihood, results)
%
% This function runs the increase-tau algorithm.

% For the notations below:
% N = no of vertices per hemisphere; 
% M = no of cortical vertices for both hemispheres;
% k = total no of parcels 

% Input
%   - params: (struct)
%     The structure containing the input arguments.
%   - label: (matrix)
%     Mx1 resultant parcellation label at the current step.
%   - full_corr_mat: (matrix)
%     The 2Nx2N premultiplied matrix computed from multiplying normalized concatenated fMRI data.
%   - lh_mask: (matrix)
%     Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the left hemisphere.
%   - rh_mask: (matrix)
%     Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the right hemisphere.
%   - avg_mesh: (struct)
%     The structure containing meshes for both left and right hemispheres.
%   - Neighborhood: (matrix)
%     An MxM matrix representing the neighborhood for both hemispheres, both interhemispheric and intrahemispheric.
%   - Smooth_cost: (matrix)
%     A kxk matrix representing smoothcost and to be fed to the GCO algorithm.
%   - tau: (matrix)
%     An 1xk array representing the current values for the tau hyperparameter.
%   - i: (integer)
%     An integer controlling the run iterations and termination criteria for GCO.
%   - likelihood: (matrix)
%     A Mxk matrix representing the current likelihood given a fixed set of labels.
%   - results: (struct)
%     Struct containing some critical intermediate information at the current optimization iteration.

% Output
%   - likelihood: (matrix)
%     The Mxk updated likelihoods after tau adaptation.
%   - results: (struct)
%     The updated results after tau adaptation.
%   - label: (struct)
%     The Mx1 updated parcellation labels after tau adaptation.

% Example
%   - [likelihood, results, label] = CBIG_hMRF_adapt_tau(params, label, full_corr_mat, lh_mask, rh_mask,...
%     avg_mesh, Neighborhood, Smooth_cost, tau, 4, likelihood, results)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
fileID = fopen(params.convergence_log_path, 'a');
fprintf(fileID, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
fprintf(fileID, 'Perform increase tau algorithm... \n');
fprintf(fileID, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
fclose(fileID);

for j = 1:params.increase_tau_iters
    fprintf('Now is the %dth iteration... \n', j);
    tau_mod = update_tau_vector(label, tau, params, avg_mesh, lh_mask, rh_mask, params.min_verts_per_cluster);
    if(tau_mod == tau)
        disp('Tau has stopped updating, no more split parcels.');
        break;
    end
    tau = tau_mod;
    
    [likelihood, results,label, tau] = CBIG_hMRF_update_labels_via_graphcut(full_corr_mat, params,...
    lh_mask, rh_mask, avg_mesh, label, Neighborhood, Smooth_cost, i, tau);
end
results.final_tau = tau;
end


function tau_mod = update_tau_vector(label, tau, params, avg_mesh, lh_mask, rh_mask, max_vertices_in_an_island)
% this function checks the existence of split parcels, then update tau for these specific parcels.
% for parcels that are split into a main component (>3 vertices) and an island (<=3 vertices), 
% we shall not increase tau for the parcel.
% for parcels that have more than 2 components (expected to be very unlikeli), and increase tau by more folds.

lh_label = label(1:sum(lh_mask));
rh_label = label(sum(lh_mask)+1:end);
lh_full_label(lh_mask) = lh_label;
rh_full_label(rh_mask) = rh_label;
[lh_ci, ~] = CBIG_hMRF_generate_components_one_hemi(avg_mesh.lh_avg_mesh, lh_full_label);
[rh_ci, ~] = CBIG_hMRF_generate_components_one_hemi(avg_mesh.rh_avg_mesh, rh_full_label);


disp('Compute number of distributed parcels (num_components - num_labels)...'); 
sp_lh = length(unique(lh_ci)) - length(unique(lh_full_label));
sp_rh = length(unique(rh_ci)) - length(unique(rh_full_label));
fprintf('lh and rh split parcels count: %d and %d \n', sp_lh, sp_rh);

disp('STEP1: Detect split parcels and their sizes on lh...');
lh_split_cluster_idx = zeros(1, params.num_cluster_per_hemi);
for which_cluster = 1:params.num_cluster_per_hemi
    if(current_cluster_need_tau_to_increase(lh_full_label, which_cluster, lh_ci, max_vertices_in_an_island))
        lh_split_cluster_idx(which_cluster) = 1;
    end
end

disp('STEP2: Detect split parcels and their sizes on rh...'); 
rh_split_cluster_idx = zeros(1, params.num_cluster_per_hemi);
for which_cluster = 1+params.num_cluster_per_hemi : params.num_cluster
    if(current_cluster_need_tau_to_increase(rh_full_label, which_cluster, rh_ci, max_vertices_in_an_island))
        rh_split_cluster_idx(which_cluster - params.num_cluster_per_hemi) = 1;
    end
end

disp('STEP3: adjust tau accordingly...'); 
lh_rh_split_cluster_idx = logical([lh_split_cluster_idx rh_split_cluster_idx]);
tau_mod = tau;
tau_mod(lh_rh_split_cluster_idx) = tau(lh_rh_split_cluster_idx) .* params.tau_increase_speed;
end

function need_increase = current_cluster_need_tau_to_increase(hemi_full_label,...
which_cluster, hemi_ci, max_vertices_in_an_island)
cluster_idx = hemi_full_label == which_cluster;
unique_comps_for_cur_cluster = unique(hemi_ci(cluster_idx));
cnt_of_comp_sizes = zeros(1,length(unique_comps_for_cur_cluster));
for comp_idx = 1:length(unique_comps_for_cur_cluster)
    cur_comp = unique_comps_for_cur_cluster(comp_idx);
    cnt_of_comp_sizes(comp_idx) = length(find(cur_comp == hemi_ci)); 
end

if(length(cnt_of_comp_sizes) > 2)
    warning('Detected a cluster of more than 2 components on left hemisphere.');
    need_increase = 1;
elseif(length(cnt_of_comp_sizes) == 2)
    if(min(cnt_of_comp_sizes) <= max_vertices_in_an_island)
        disp('Current split parcel contains a main component and an island, skipping now.');
        disp('The island will be absorbed by the function CBIG_absorb_tiny_clusters subsequently');
        need_increase = 0;
    else
        fprintf('Sizes of components of the split parcels: %d  %d \n', cnt_of_comp_sizes(1), cnt_of_comp_sizes(2));
        need_increase = 1;
    end
else
    need_increase = 0;
end
end