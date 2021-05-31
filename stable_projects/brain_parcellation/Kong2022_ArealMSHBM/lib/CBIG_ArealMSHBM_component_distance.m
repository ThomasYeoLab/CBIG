function eucli_dist = CBIG_ArealMSHBM_component_distance(lh_labels, rh_labels, mesh, num_parcel)
    
% eucli_dist = CBIG_ArealMSHBM_component_distance(lh_labels, rh_labels, mesh, num_parcel)
%
% This script will compute the Euclidean distances between distributed components for each 
% parcel of a given parcellation <lh_labels> <rh_labels>. The output <eucli_dist> is defined
% as the maximum shortest distance among components. Parcels with a single component (i.e. 
% contiguous parcels) will be denoted as 0 in the output vector <eucli_dist>. For example, if
% parcel 1 has 3 components A,B,C. The shortest distance between A and B is 3, the shortest 
% distance between A and C is 2, the shortest distance between B and C is 4. Then eucli_dist(1)
% = 4.
%
% Input:
%   - lh_labels, rh_labels: (Nx1 vector)
%     The given parcellation labels. Medial wall should be denoted as 0.
%
%   - mesh: (string)
%     The mesh name of the parcellation surface space. For example, 'fsaverage6', 'fs_LR_32k'.
%
%   - num_parcel: (scalar)
%     The total number of parcels (excluding medial wall). For example, 400.
%
% Output:
%   - eucli_dist: (1xnum_parcel vector)
%     The maximum shortest distances among components for each parcel of a given parcellation.
%
% Example:
%   eucli_dist = CBIG_ArealMSHBM_component_distance(lh_labels, rh_labels, 'fs_LR_32k', 400);
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib'));

%% Find distributed parcels
parcel_components = CBIG_ArealMSHBM_compute_components_general(lh_labels, rh_labels, mesh, num_parcel);
parcels = find(parcel_components == 1);
parcels_dist = find(parcel_components > 1);
eucli_dist = zeros(size(parcel_components));

%% If there exists any distributed parcel
if(~isempty(parcels_dist))
    for i = parcels
        lh_labels(lh_labels == i) = 0;
        rh_labels(rh_labels == i) = 0;
    end
    rh_labels(rh_labels~=0) = rh_labels(rh_labels~=0) - num_parcel/2;
    if(strcmp(mesh, 'fs_LR_32k'))
        lh_avg_mesh=CBIG_read_fslr_surface('lh','fs_LR_32k','sphere');
        rh_avg_mesh=CBIG_read_fslr_surface('rh','fs_LR_32k','sphere');
    else
        lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'sphere','cortex');
        rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'sphere', 'cortex');
    end
    % Find distributed componnets for distributed parcels
    [lh_ci,lh_sizes]=CBIG_ArealMSHBM_find_components(lh_avg_mesh, lh_labels);
    [rh_ci,lh_sizes]=CBIG_ArealMSHBM_find_components(rh_avg_mesh, rh_labels);
    lh_labels_tmp = lh_labels;
    rh_labels_tmp = rh_labels;

    lh_labels_tmp(lh_labels == 0) = num_parcel + 1;
    rh_labels_tmp(rh_labels == 0) = num_parcel + 1;

    % Find boundary vertices for all components
    [lh_up_labels,rh_up_labels]=CBIG_ArealMSHBM_BuildTwoVertThickBoundary(lh_avg_mesh,...
        rh_avg_mesh,lh_labels_tmp,rh_labels_tmp);
    lh_labels_bound = lh_labels.*(lh_up_labels==0);
    rh_labels_bound = rh_labels.*(rh_up_labels==0);

    % Compute Euclidean distances among boudary vertices
    lh_verts_xyz = lh_avg_mesh.vertices(:, lh_labels_bound~=0);
    rh_verts_xyz = rh_avg_mesh.vertices(:, rh_labels_bound~=0);
    lh_pdist = squareform(pdist(lh_verts_xyz'));
    rh_pdist = squareform(pdist(rh_verts_xyz'));
    lh_labels_bound(lh_labels_bound == 0) = [];
    rh_labels_bound(rh_labels_bound == 0) = [];
    lh_bound_mat = zeros(length(lh_labels_bound), length(parcel_components)./2);
    rh_bound_mat = zeros(length(rh_labels_bound), length(parcel_components)./2);
    lh_bound_mat(sub2ind(size(lh_bound_mat),[1:length(lh_labels_bound)]',lh_labels_bound)) = 1;
    lh_bound_mask = lh_bound_mat*lh_bound_mat';
    rh_bound_mat(sub2ind(size(rh_bound_mat),[1:length(rh_labels_bound)]',rh_labels_bound)) = 1;
    rh_bound_mask = rh_bound_mat*rh_bound_mat';

    lh_ci(lh_labels == 0) = 0;
    lh_comp_bound = lh_ci.*(lh_up_labels==0);
    lh_comp_bound(lh_comp_bound == 0) = [];
    lh_comp_bound_mat = zeros(length(lh_comp_bound), max(lh_comp_bound));
    lh_comp_bound_mat(sub2ind(size(lh_comp_bound_mat),[1:length(lh_comp_bound)]',lh_comp_bound)) = 1;
    lh_comp_bound_mask = lh_comp_bound_mat*lh_comp_bound_mat';
    lh_pdist_mask = lh_pdist.*lh_bound_mask.*(lh_comp_bound_mask==0);

    rh_ci(rh_labels == 0) = 0;
    rh_comp_bound = rh_ci.*(rh_up_labels==0);
    rh_comp_bound(rh_comp_bound == 0) = [];
    rh_comp_bound_mat = zeros(length(rh_comp_bound), max(rh_comp_bound));
    rh_comp_bound_mat(sub2ind(size(rh_comp_bound_mat),[1:length(rh_comp_bound)]',rh_comp_bound)) = 1;
    rh_comp_bound_mask = rh_comp_bound_mat*rh_comp_bound_mat';
    rh_pdist_mask = rh_pdist.*rh_bound_mask.*(rh_comp_bound_mask==0);

    % Single component parcels will have 0 distance
    lh_pdist_mask(lh_pdist_mask==0) = Inf;
    lh_min_bound_pdist = min(lh_pdist_mask);
    lh_min_bound_pdist(isinf(lh_min_bound_pdist)) = 0;
    lh_comp_dist = lh_min_bound_pdist'.*lh_comp_bound_mat;
    lh_comp_dist(lh_comp_dist==0) = Inf;
    lh_comp_min = min(lh_comp_dist);
    lh_comp_min(isinf(lh_comp_min)) = 0;
    
    % Map components to parcels
    if(~isempty(lh_comp_bound_mat))
        lh_map_mat = lh_bound_mat'*lh_comp_bound_mat ~= 0;
        lh_max_min_dist = max(lh_comp_min.*lh_map_mat,[],2);
    else
        lh_max_min_dist = zeros(size(parcel_components,1),size(parcel_components,2)/2)';
    end
    
    % Single component parcels will have 0 distance
    rh_pdist_mask(rh_pdist_mask==0) = Inf;
    rh_min_bound_pdist = min(rh_pdist_mask);
    rh_min_bound_pdist(isinf(rh_min_bound_pdist)) = 0;
    rh_comp_dist = rh_min_bound_pdist'.*rh_comp_bound_mat;
    rh_comp_dist(rh_comp_dist==0) = Inf;
    rh_comp_min = min(rh_comp_dist);
    rh_comp_min(isinf(rh_comp_min)) = 0;

    % Map components to parcel
    if(~isempty(rh_comp_bound_mat))
        rh_map_mat = rh_bound_mat'*rh_comp_bound_mat ~= 0;
        rh_max_min_dist = max(rh_comp_min.*rh_map_mat,[],2);
    else
        rh_max_min_dist = zeros(size(parcel_components,1),size(parcel_components,2)/2)';
    end

    eucli_dist = [lh_max_min_dist;rh_max_min_dist];
end
eucli_dist = eucli_dist';

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib'));

end
