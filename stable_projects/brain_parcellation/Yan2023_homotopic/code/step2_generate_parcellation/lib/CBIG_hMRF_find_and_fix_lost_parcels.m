function [label, tau] = CBIG_hMRF_find_and_fix_lost_parcels(label, avg_mesh, lh_mask, rh_mask,...
     full_corr_mat, tau, params)
% [label, tau] = CBIG_hMRF_find_and_fix_lost_parcels(label, avg_mesh, lh_mask, rh_mask, full_corr_mat, tau, params)
%
% This function processes lost parcels given a set of parcellation labels.
% Lost parcels are defined as those that are missing on one or both hemispheres.
% First the function merges those singular parcels (exist on one hemisphere only.
% Then the function look for lost parcel pairs (non existent on either hemispheres.)
% Afterwards it generate a list of existing parcels pairs that are suitable for being split into 2 parcel pairs.
% It will then assign the lost parcel labels into newly generated parcel pairs.

% For the notations below:
% N = no of vertices per hemisphere; 
% M = no of cortical vertices for both hemispheres;
% k = total no of parcels 

% Input
%   - label: (matrix)
%     Resultant Mx1 parcellation label at the current step.
%   - avg_mesh: (struct)
%     The structure containing meshes for both left and right hemispheres.
%   - lh_mask: (matrix)
%     A Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the left hemisphere.
%   - rh_mask: (matrix)
%     A Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the right hemisphere.
%   - full_corr_mat: (matrix)
%     The 2Nx2N premultiplied matrix computed from multiplying normalized concatenated fMRI data.
%   - tau: (matrix)
%     An 1xk array representing the current values for the tau hyperparameter.
%   - params: (struct)
%     The structure containing the input arguments.
% Output
%   - label: (matrix)
%     Resultant Mx1 parcellation label at the current step.
%   - tau: (matrix)
%     An 1xk array representing the current values for the tau hyperparameter. Might differ from initial input tau.
%
% Example
%   - [label, tau] = CBIG_hMRF_find_and_fix_lost_parcels(label, avg_mesh, lh_mask, rh_mask, full_corr_mat, tau, params)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

disp('##############################');
disp('# Processing lost parcels... #');
disp('##############################');

lh_avg_mesh = avg_mesh.lh_avg_mesh;
rh_avg_mesh = avg_mesh.rh_avg_mesh;
num_cortical_verts = length(lh_mask);
num_clusters_per_hemi = params.num_cluster_per_hemi;

[lh_full_label, rh_full_label] = CBIG_hMRF_get_left_right_overlapping_labels(lh_mask, rh_mask,...
    label, num_clusters_per_hemi);

disp('Checking if there are singular parcels on both hemispheres...');
[exist_on_lh_only, exist_on_rh_only] = profile_lost_parcels(lh_full_label, rh_full_label,...
    lh_mask, rh_mask, num_clusters_per_hemi);

if(~isempty(exist_on_lh_only))
    fprintf('Parcels that exist on left hemisphere only: %d \n', exist_on_lh_only);
    lh_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(full_corr_mat(1:num_cortical_verts,...
    1:num_cortical_verts), exist_on_lh_only, lh_full_label, lh_avg_mesh, lh_mask);
end
if(~isempty(exist_on_rh_only))
    fprintf('Parcels that exist on right hemisphere only: %d \n', exist_on_rh_only);
    rh_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(full_corr_mat(num_cortical_verts+1:end,...
        num_cortical_verts+1:end), exist_on_rh_only, rh_full_label, rh_avg_mesh, rh_mask);
end

% check if there are still singular parcels on either hemisphere
[exist_on_lh_only, exist_on_rh_only] = profile_lost_parcels(lh_full_label,rh_full_label, lh_mask, rh_mask,...
num_clusters_per_hemi);
assert(isempty(exist_on_lh_only), 'lh still have singular parcels.');
assert(isempty(exist_on_rh_only), 'rh still have singular parcels.');

disp('Checking if there are parcels missing on both hemispheres...');
[~, ~, lh_rh_both_missing, ~, ~] = profile_lost_parcels(lh_full_label,rh_full_label, lh_mask, rh_mask,...
num_clusters_per_hemi);

if(~isempty(lh_rh_both_missing))
    fprintf('Parcels that are missing on both hemispheres: %d \n', lh_rh_both_missing);
    [lh_full_label, rh_full_label] = assign_lost_parcel_both_hemi(full_corr_mat, lh_full_label,...
    rh_full_label, lh_mask, rh_mask, avg_mesh, num_clusters_per_hemi, lh_rh_both_missing, params);
end


% check if there are still any singular parcel on either hemisphere, or lost on both hemispheres
[exist_on_lh_only, exist_on_rh_only, lh_rh_both_missing] = profile_lost_parcels(lh_full_label,rh_full_label,...
lh_mask, rh_mask, num_clusters_per_hemi);
assert(isempty(exist_on_lh_only), 'lh still have singular parcels.');
assert(isempty(exist_on_rh_only), 'rh still have singular parcels.');
assert(isempty(lh_rh_both_missing), 'There are parcels lost on both hemispheres that remain unassigned.');

% for the newly assigned empty clusters, need to restore tau to the original value
tau(exist_on_lh_only) = params.initial_tau;
tau(exist_on_rh_only) = params.initial_tau;
tau(lh_rh_both_missing) = params.initial_tau;

rh_full_label = rh_full_label + double(num_clusters_per_hemi);
label = [lh_full_label(lh_mask) rh_full_label(rh_mask)];
end

function [on_lh_only, on_rh_only, lh_rh_both_missing, lh_target_parcel_sizes, rh_target_parcel_sizes] =...
    profile_lost_parcels(lh_full_label, rh_full_label,lh_mask, rh_mask, num_of_clusters)
% this function identifies the lost parcels and count the sizes for the singular parcels

lh_label = lh_full_label(lh_mask);
rh_label = rh_full_label(rh_mask);
lh_label_unique = unique(lh_label);
rh_label_unique = unique(rh_label);

lh_missing = [];
rh_missing = [];

if(length(lh_label_unique) ~= num_of_clusters)
    lh_missing = setdiff(1: num_of_clusters, lh_label_unique);
end

if(length(rh_label_unique) ~= num_of_clusters)
    rh_missing = setdiff(1: num_of_clusters, rh_label_unique);
end

% case I: lost on both hemi
lh_rh_both_missing = intersect(lh_missing, rh_missing);
% case II: lost on one hemi
on_lh_only = setdiff(lh_label_unique, rh_label_unique);
on_rh_only = setdiff(rh_label_unique, lh_label_unique);

% count sizes
[lh_parcel_size_count, ~] = histcounts(lh_label, 1:num_of_clusters+1);
[rh_parcel_size_count, ~] = histcounts(rh_label, 1:num_of_clusters+1);

% for those that exist only on lh
lh_target_parcel_sizes = lh_parcel_size_count(on_lh_only);

% for those that exist only on rh
rh_target_parcel_sizes = rh_parcel_size_count(on_rh_only);

end


function [lh_full_label_mod, rh_full_label_mod] = assign_lost_parcel_both_hemi(full_lh_rh_corr,...
    lh_full_label, rh_full_label, lh_mask, rh_mask, avg_mesh, num_clusters_per_hemi, lost_parcel_list, params)
% this function reassign a list of lost parcels for both hemispheres concurrently.

if(size(lh_mask,1) < size(lh_mask,2))
lh_mask = lh_mask';
end

if(size(rh_mask,1) < size(rh_mask,2))
    rh_mask = rh_mask';
end

full_mask = [lh_mask; rh_mask];

min_verts_candidate_parcel = floor(sum(full_mask)/(num_clusters_per_hemi*4));

disp('STEP0: Prep work...');
[lh_ci, ~] = CBIG_hMRF_generate_components_one_hemi(avg_mesh.lh_avg_mesh, lh_full_label);
[rh_ci, ~] = CBIG_hMRF_generate_components_one_hemi(avg_mesh.rh_avg_mesh, rh_full_label);

disp('STEP1: Find candidate cluster list...');
% homogeneity will be averaged for pairs of components on lh/rh 
% clusters that have parcel loss / split parcels will be filtered out from the list
candidate_clusters = find_candidate_clusters(lh_full_label, rh_full_label,...
    num_clusters_per_hemi, lh_ci, rh_ci, lh_mask, rh_mask, min_verts_candidate_parcel);

disp('STEP2: Order the k clusters by decreasing homogeneity...');
candidate_cluster_ascending_homo = order_candidates_by_homogeneity(full_lh_rh_corr,...
candidate_clusters, lh_full_label, rh_full_label, lh_ci, rh_ci);

disp('STEP3: Split the candidate cluster, reassign for lost parcel...');
[lh_full_label_mod, rh_full_label_mod] = divide_parcel_and_reassign(lost_parcel_list,...
candidate_cluster_ascending_homo, lh_full_label, rh_full_label, full_lh_rh_corr, avg_mesh, params);
end

function [lh_full_label, rh_full_label] = divide_parcel_and_reassign(lost_parcel_list,...
candidate_cluster_ascending_homo, lh_full_label, rh_full_label, full_lh_rh_corr, avg_mesh, params)
% while there still remains lost parcels to re-assign, we would keep dividing existing homotopic pairs
% into two smaller parcel pairs and try to reassign one of these to the lost parcel label

candidate_idx = 1;
k_means_rand_inits = 3;
for i = 1:length(lost_parcel_list)
    lost_parcel_idx = lost_parcel_list(i);
    % keep popping from candidate_cluster_ascending_homo 
    % in case the current one kept divide into disconnected parcels
    while(true)
        candidate_parcel_idx = candidate_cluster_ascending_homo(candidate_idx);
        candidate_idx = candidate_idx + 1;
        if(candidate_idx > length(candidate_cluster_ascending_homo))
            warning('There are not enough parcels to be split and reassigned. Please adjust the hyperparameters.');
            break;
        end

        [lh_candidate_pair_list] = divide_single_parcel(lh_full_label,...
            full_lh_rh_corr(1:length(lh_full_label),:), avg_mesh.lh_avg_mesh, candidate_parcel_idx,...
            k_means_rand_inits);
        [rh_candidate_pair_list] = divide_single_parcel(rh_full_label,...
            full_lh_rh_corr(length(lh_full_label)+1:end,:), avg_mesh.rh_avg_mesh, candidate_parcel_idx,...
            k_means_rand_inits);
        
        if(~isempty(lh_candidate_pair_list) && ~isempty(rh_candidate_pair_list))
            break;
        end
    end
    
    [lh_full_label, rh_full_label] = pick_best_pair_from_candidates_and_reassign(lh_candidate_pair_list,...
        rh_candidate_pair_list, lh_full_label, rh_full_label, lost_parcel_idx, params);
end

end

function candidate_cluster_ascending_homo = order_candidates_by_homogeneity(full_lh_rh_corr,...
candidate_clusters, lh_full_label, rh_full_label, lh_ci, rh_ci)
% order candidate parcels by ascending homogeneity
num_cortical_verts = length(lh_full_label);
lh_full_corr = full_lh_rh_corr(1:num_cortical_verts, 1:num_cortical_verts);
rh_full_corr = full_lh_rh_corr(num_cortical_verts+1:end, num_cortical_verts+1:end);

num_candidates = length(candidate_clusters);
lh_homo_mat = zeros(num_candidates,2);
rh_homo_mat = zeros(num_candidates,2);

for i = 1: num_candidates
    lh_ci_k = lh_ci(lh_full_label == (candidate_clusters(i)));
    if(length(unique(lh_ci_k)) == 1)
        temp = lh_full_corr(lh_full_label==candidate_clusters((i)), lh_full_label==candidate_clusters(i));
        lh_homo_mat(i,:) = [candidate_clusters(i), (sum(sum(temp))-size(temp,1))/(size(temp,1)*(size(temp,1)-1))];
    end
    rh_ci_k = rh_ci(rh_full_label == candidate_clusters(i));
    if(length(unique(rh_ci_k)) == 1)
        temp = rh_full_corr(rh_full_label==candidate_clusters(i), rh_full_label==candidate_clusters(i));
        rh_homo_mat(i,:) = [candidate_clusters(i), (sum(sum(temp))-size(temp,1))/(size(temp,1)*(size(temp,1)-1))];
    end
end

lh_homo_mat = lh_homo_mat(any(lh_homo_mat,2), :);
rh_homo_mat = rh_homo_mat(any(rh_homo_mat,2), :);

lhrh_homo_mat(:,1) = lh_homo_mat(:,1);
lhrh_homo_mat(:,2) = (lh_homo_mat(:,2) + rh_homo_mat(:,2))/2;
lhrh_homo_mat = sortrows(lhrh_homo_mat, 2);
candidate_cluster_ascending_homo = lhrh_homo_mat(:,1);
end

function candidate_clusters = find_candidate_clusters(lh_full_label, rh_full_label,...
num_clusters_per_hemi, lh_ci, rh_ci, lh_mask, rh_mask, min_verts_candidate_parcel)
% find candidate symmetric parcel pairs to split into two homotopic pairs

% find which cluster appear symmetrically on both hemispheres
lh_rh_symmetric_clusters_idx = find_symmetric_cluster_pairs(lh_full_label, rh_full_label, num_clusters_per_hemi);

% find which clusters contain split parcels
lh_rh_split_cluster_idx = detect_split_parcels(lh_full_label, rh_full_label, lh_ci, rh_ci, num_clusters_per_hemi);

% count the parcel sizes, find the clusters that are big enough to be split and reassigned
lh_label = lh_full_label(lh_mask);
rh_label = rh_full_label(rh_mask);
lh_rh_big_cluster_idx = find_big_clusters(lh_label, rh_label, num_clusters_per_hemi, min_verts_candidate_parcel);

% find intersection of the three vectors
final_candidate_cluster_idx = lh_rh_symmetric_clusters_idx' & ~lh_rh_split_cluster_idx & lh_rh_big_cluster_idx;
candidate_clusters = find(final_candidate_cluster_idx == 1);

end

function big_cluster_idx = find_big_clusters(lh_label, rh_label, num_clusters_per_hemi, min_verts_candidate_parcel)
% find big parcels by ranking parcels by number of enclosed vertices

[lh_hist_counts,~] = histcounts(lh_label,1:num_clusters_per_hemi+1);
lh_big_cluster_idx = lh_hist_counts >= min_verts_candidate_parcel;
[rh_hist_counts,~] = histcounts(rh_label,1:num_clusters_per_hemi+1);
rh_big_cluster_idx = rh_hist_counts >= min_verts_candidate_parcel;
big_cluster_idx = lh_big_cluster_idx | rh_big_cluster_idx;

end

function lh_rh_symmetric_clusters_idx = find_symmetric_cluster_pairs(lh_full_label, rh_full_label,...
    num_clusters_per_hemi)
% find which cluster appear symmetrically on both hemispheres

% intersect lh and rh labels; filter out parcels that have lost counterparts
lh_rh_symmetric_clusters = nonzeros(intersect(lh_full_label, rh_full_label));
lh_rh_symmetric_clusters_idx = zeros(num_clusters_per_hemi, 1);
lh_rh_symmetric_clusters_idx(lh_rh_symmetric_clusters) = 1;

end

function [lh_full_label, rh_full_label] = pick_best_pair_from_candidates_and_reassign(lh_candidate_pair_list,...
    rh_candidate_pair_list, lh_full_label, rh_full_label, lost_parcel_idx, params)
% pick the most suitable candidate parcel pairs and reassign parcel labels

best_mean_overlap = -1e10;
for lh_idx = 1:length(lh_candidate_pair_list)
    for rh_idx = 1:length(rh_candidate_pair_list)
        lh_candidate_pair = lh_candidate_pair_list{lh_idx};
        rh_candidate_pair = rh_candidate_pair_list{rh_idx};

        [cur_lh_rh_mean_overlap, is_match_A_to_A] = match_candidate_pairs_by_max_overlap(lh_candidate_pair,...
            rh_candidate_pair, params);
        
        if(cur_lh_rh_mean_overlap > best_mean_overlap)
            best_mean_overlap = cur_lh_rh_mean_overlap;
            lh_idx_best = lh_idx;
            rh_idx_best = rh_idx;
            is_match_A_to_A_best = is_match_A_to_A;
        end
    end 
end  

best_lh_candidate = lh_candidate_pair_list{lh_idx_best};
best_rh_candidate = rh_candidate_pair_list{rh_idx_best};

[lh_full_label, rh_full_label] = reassign_labels_for_best_candidate(best_lh_candidate,...
    best_rh_candidate, is_match_A_to_A_best, lh_full_label, rh_full_label, lost_parcel_idx);

end

function [lh_full_label, rh_full_label] = reassign_labels_for_best_candidate(best_lh_candidate,...
best_rh_candidate, is_match_A_to_A_best, lh_full_label, rh_full_label, lost_parcel_idx)
% perform label reassignment

lh_A_verts = best_lh_candidate.parcel_half_A_verts_idx;

rh_A_verts = best_rh_candidate.parcel_half_A_verts_idx;
rh_B_verts = best_rh_candidate.parcel_half_B_verts_idx;

if(is_match_A_to_A_best)
    lh_full_label(lh_A_verts) = lost_parcel_idx;
    rh_full_label(rh_A_verts) = lost_parcel_idx;
else
    lh_full_label(lh_A_verts) = lost_parcel_idx;
    rh_full_label(rh_B_verts) = lost_parcel_idx;
end

end

function [cur_lh_rh_mean_overlap,is_match_A_to_A] = match_candidate_pairs_by_max_overlap(lh_candidate_pair,...
rh_candidate_pair, params)
% lh and rh parcel both have parts A and B; match across lh and rh by comparing % of overlap

lh_A_verts = lh_candidate_pair.parcel_half_A_verts_idx;
rh_A_verts = rh_candidate_pair.parcel_half_A_verts_idx;
rh_B_verts = rh_candidate_pair.parcel_half_B_verts_idx;

perc_overlapped_A2A = count_percentage_overlapped(lh_A_verts, rh_A_verts, params);
perc_overlapped_A2B = count_percentage_overlapped(lh_A_verts, rh_B_verts, params);

if(perc_overlapped_A2A > perc_overlapped_A2B)
    is_match_A_to_A = 1;
    cur_lh_rh_mean_overlap = perc_overlapped_A2A;
else
    is_match_A_to_A = 0;
    cur_lh_rh_mean_overlap = perc_overlapped_A2B;
end
end

function perc_overlapped = count_percentage_overlapped(lh_verts, rh_verts, params)
% compute parcel cross-hemispheric overlap by checking percentage of vertices that are neighbors in the 
% cross-hemispheric neighborhood

load(params.lh_vert2rh_vert, 'lh_one2rh_one');
count_matched_lhrh = 0;

for idx = 1:length(lh_verts)
    cur_lh_vert = lh_verts(idx);
    cur_lh_vert_matched_to_rh = lh_one2rh_one(cur_lh_vert);
    if(ismember(cur_lh_vert_matched_to_rh, rh_verts))
        count_matched_lhrh = count_matched_lhrh + 1;
    end
end
perc_overlapped = count_matched_lhrh / (length(lh_verts) + length(rh_verts));
end

function lh_rh_split_cluster_idx = detect_split_parcels(lh_full_label, rh_full_label, lh_ci, rh_ci, num_cluster)
% this function detect split parcels on both hemispheres

disp('Compute number of distributed parcels (num_components - num_labels)...'); 
sp_lh = length(unique(lh_ci)) - length(unique(lh_full_label));
sp_rh = length(unique(rh_ci)) - length(unique(rh_full_label));
fprintf('lh and rh split parcels count: %d and %d \n', sp_lh, sp_rh);

disp('STEP1: Detect split parcels on lh...');
lh_split_cluster_idx = zeros(1, num_cluster);
for which_cluster = 1:num_cluster
    cluster_verts_idx = lh_full_label == which_cluster;
    unique_comp_for_cur_cluster = unique(lh_ci(cluster_verts_idx));
    if(length(unique_comp_for_cur_cluster) > 1)
        lh_split_cluster_idx(which_cluster) = 1;
    end
end    

disp('STEP2: Detect split parcels on rh...');
rh_split_cluster_idx = zeros(1, num_cluster);
for which_cluster = 1:num_cluster
    cluster_verts_idx = rh_full_label == which_cluster;
    unique_comp_for_cur_cluster = unique(rh_ci(cluster_verts_idx));
    if(length(unique_comp_for_cur_cluster) > 1)
        rh_split_cluster_idx(which_cluster) = 1;
    end
end

lh_rh_split_cluster_idx = (lh_split_cluster_idx | rh_split_cluster_idx);
end

function [candidate_pairs] = divide_single_parcel(hemi_full_label, hemi_full_corr,...
    hemi_mesh, candidate_parcel_idx, rand_inits)
% split parcel using kmeans into 2 halves

parcel_verts_idx = find(hemi_full_label == candidate_parcel_idx);
% for each random initialization, we need to keep record of the following properties
candidate_pairs = {};

for iteration = 1:rand_inits
    rng(iteration, 'philox');
    [divided_parcel_label_assignment, ~, sumd] = kmeans(hemi_full_corr(parcel_verts_idx, :), 2);

    pair.parcel_original_idx = candidate_parcel_idx;
    pair.parcel_half_A_verts_idx = parcel_verts_idx(divided_parcel_label_assignment == 1);
    pair.parcel_half_B_verts_idx = parcel_verts_idx(divided_parcel_label_assignment == 2);
    pair.sumd = sumd;
    % assume we have already divided the parcel, need to check if the
    % division introduces new disconnected parcels
    hemi_full_label_region_modified = zeros(size(hemi_full_label));
    hemi_full_label_region_modified(pair.parcel_half_A_verts_idx) = 1;
    hemi_full_label_region_modified(pair.parcel_half_B_verts_idx) = 2;
    
    if(~region_contain_disconnected_parcels(hemi_mesh, hemi_full_label_region_modified))
        candidate_pairs{end+1} = pair;
    end
end

end

function is_split = region_contain_disconnected_parcels(hemi_mesh, hemi_full_label_region_modified)
% check if the given label contains disconnected parcels

[ci,~] = CBIG_hMRF_generate_components_one_hemi(hemi_mesh, hemi_full_label_region_modified);
if(length(unique(ci)) > length(unique(hemi_full_label_region_modified)))
    is_split = true;
    return
else
    is_split = false;
    return
end
end