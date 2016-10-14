function CBIG_ClassifyBrainVolumeBasedOnSurf2VolProfile(input_profile, mask_file, input_surf_clusters, output_prefix, lh_input_surf_profile, rh_input_surf_profile)

% CBIG_ClassifyBrainVolumeBasedOnSurf2VolProfile(input_profile, mask_file, input_surf_clusters, output_prefix, lh_input_surf_profile, rh_input_surf_profile)
% 
% Input arguments:
%   - input_profile         : surf2vol correlation profile, which could be
%                             created by CBIG_ComputeCorrelationProfileSurf2Vol*.m
%   - mask_file             : a user specified volume mask
%   - input_surf_clusters   : surface clustering result (.mat)
%   - output_prefix         : output file name without extensions
%   - lh_input_surf_profile : surface correlation profile on left hemisphere
%   - rh_input_surf_profile : surface correlation profile on right hemisphere
%
% Classify volumetric voxels based on surface to volume correlation
% profiles and pre-computed surface parcellation result. The volumetric
% parcellation, its log probability, as well as silhouette values are saved
% out.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% read mask
disp(['Reading mask: ' mask_file]);
mask = MRIread(mask_file);
mask_index = find(mask.vol == 1);

% read profile
disp(['Reading profile: ' input_profile]);
load(input_profile);
series = surf2vol_correlation_profile;
clear surf2vol_correlation_profile;

% von_mises pre-filtering
series = series - repmat(mean(series, 2), 1, size(series, 2));
series = series ./ repmat(sqrt(sum(series.^2, 2)),1,size(series, 2));

% read surface clusters
disp(['Reading surface clusters: ' input_surf_clusters]);
surf_clusters = load(input_surf_clusters);
lh_labels = surf_clusters.lh_labels;
rh_labels = surf_clusters.rh_labels;

% classify brain volume
num_clusters = size(surf_clusters.mtc, 1);
distctrs = surf_clusters.lambda*series*(surf_clusters.mtc');
maxdists = max(distctrs,[],2);

rr = exp(distctrs - maxdists*ones(1, num_clusters));
sumrr = sum(rr,2);
prob = rr ./ sumrr(:,ones(num_clusters,1));

[maxr, assign] = max(prob,[],2);

% write out results
output_file = [output_prefix '.nii.gz'];
disp(['Writing parcellation: ' output_file]);
output = mask;
output.vol(mask_index) = assign;
MRIwrite(output, output_file);

output_file = [output_prefix '.logprob.nii.gz'];
disp(['Writing parcellation: ' output_file]);
size_input = size(mask.vol);
size_input(4) = num_clusters;
output.vol = zeros(size_input);
frame = zeros(size_input(1:3));
logprob = distctrs - maxdists*ones(1, num_clusters);
for i = 1:num_clusters
    frame(mask_index) = logprob(:, i);
    output.vol(:, :, :, i) = frame;
end
disp(output_file);
MRIwrite(output, output_file);

% Compute silhouette
disp('Computing Silhouette');
if(num_clusters > 1)
    tic
    disp(['Reading lh surf profiles: ' lh_input_surf_profile]);
    lh_surf_profile = MRIread(lh_input_surf_profile);
    disp(['Reading rh surf profiles: ' rh_input_surf_profile]);
    rh_surf_profile = MRIread(rh_input_surf_profile);
    lh_surf_profile = reshape(lh_surf_profile.vol, size(lh_surf_profile.vol, 1)*size(lh_surf_profile.vol, 2)*size(lh_surf_profile.vol, 3), size(lh_surf_profile.vol, 4));
    rh_surf_profile = reshape(rh_surf_profile.vol, size(rh_surf_profile.vol, 1)*size(rh_surf_profile.vol, 2)*size(rh_surf_profile.vol, 3), size(rh_surf_profile.vol, 4));

    surf_profile = [lh_surf_profile(lh_labels~=0, :); rh_surf_profile(rh_labels~=0, :)]; % voxels x ROI
    labels = [lh_labels(lh_labels~=0); rh_labels(rh_labels~=0)]; % voxels x 1

    surf_profile = surf_profile - repmat(mean(surf_profile, 2), 1, size(surf_profile, 2));
    surf_profile = surf_profile ./ repmat(sqrt(sum(surf_profile.^2, 2)), 1, size(surf_profile, 2));

    s = zeros(size(series, 1), 1);
    for i = 1:length(mask_index)

        corr_val = sum(repmat(series(i, :), size(surf_profile, 1), 1) .* surf_profile, 2);
        cluster_dist = zeros(max(assign), 1)+inf;
        for j = 1:max(assign)
            if(~isempty(assign == j))
                cluster_dist(j) = 1 - mean(corr_val(labels == j));
            end
        end
        within_clust_dist = cluster_dist(assign(i));
        cluster_dist(assign(i)) = inf;
        closest_clust_dist = min(cluster_dist);
        s(i) = (closest_clust_dist - within_clust_dist)/max(within_clust_dist, closest_clust_dist);
    end
    toc
    output_file = [output_prefix '.cort_silhouette.nii.gz'];
    output = mask;
    output.vol(mask_index) = s;
    disp(output_file);
    MRIwrite(output, output_file);

    tic; s = silhouette(series, assign, 'correlation'); toc
    output_file = [output_prefix '.subcort_silhouette.nii.gz'];
    output = mask;
    output.vol(mask_index) = s;
    disp(output_file);
    MRIwrite(output, output_file);

else
    s = zeros(size(series, 1), 1);

    output_file = [output_prefix '.cort_silhouette.nii.gz'];
    output = mask;
    output.vol(mask_index) = s;
    disp(output_file);
    MRIwrite(output, output_file);

    output_file = [output_prefix '.subcort_silhouette.nii.gz'];
    disp(output_file);
    MRIwrite(output, output_file);
end




