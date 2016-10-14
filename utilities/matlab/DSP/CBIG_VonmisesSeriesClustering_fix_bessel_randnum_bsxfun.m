function CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, num_tries, normalize, max_iter)

% CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, num_tries, normalize, max_iter)
%
% Von Mises-Fisher clustering on surface data.
% 
% Input arguments:
%     - mesh_name    : e.g. 'fsaverage5'
%     - mask         : e.g. 'cortex'
%     - num_clusters : number of clusters
%     - output_file  : output file name
%     - profile1     : group average profile on left hemisphere
%     - profile2     : group average profile on right hemisphere
%     - num_smooth   : how many times that smoothing is performed
%     - num_tries    : number of difference random initialization
%     - normailize   : 0 or 1, whether z-normalization is performed across vertices
%     - max_iter     : maximum number of iterations for one random initialization
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md




if(~exist('max_iter', 'var'))
   max_iter = 100; 
else
   if(ischar(max_iter))
       max_iter = str2num(max_iter);
   end
end


% % We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
% rand('twister',5489)

if(ischar(normalize))
   normalize = str2num(normalize); 
end

if(ischar(num_clusters))
   num_clusters = str2num(num_clusters); 
end

if(ischar(num_tries))
   num_tries = str2num(num_tries); 
end

if(ischar(num_smooth))
   num_smooth = str2num(num_smooth); 
end

% read mask
lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_name, 'inflated', mask);
l1 = find(lh_avg_mesh.MARS_label == 2); lh_num_verts = size(lh_avg_mesh.vertices, 2);
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_name, 'inflated', mask);
l2 = find(rh_avg_mesh.MARS_label == 2); rh_num_verts = size(rh_avg_mesh.vertices, 2);

% read data (voxels x N subjects)
tmp1 = MRIread(profile1);
tmp2 = MRIread(profile2);
series = [reshape(tmp1.vol, lh_num_verts, size(tmp1.vol, 4)); reshape(tmp2.vol, lh_num_verts, size(tmp2.vol, 4))]; 


% We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
% rand('twister',5489) was moved to after MRIread is because MRIread calls/apps/arch/Linux_x86_64/freesurfer/4.5.0/matlab/load_nifti.m, which calls rand('state', sum(100*clock));
rand('twister',5489)

% smooth
series(1:end/2, :)     = transpose(MARS_AverageData(lh_avg_mesh, transpose(series(1:end/2, :)), 0, num_smooth));
series(end/2+1:end, :) = transpose(MARS_AverageData(rh_avg_mesh, transpose(series(end/2+1:end, :)), 0, num_smooth));

% extract mask voxels series
series = series([l1 l2+length(lh_avg_mesh.MARS_label)], :);

% remove voxels that are completely not correlated with any rois. 
non_zero_corr_index = (sum(series, 2) ~= 0);
series = series(non_zero_corr_index, :);

% znormalize (series assumed to be voxels x subjects or voxels x profile)
if(normalize)
    mean_series = nanmean(series, 1);
    std_series = nanstd(series, 1, 1);
    series = bsxfun(@minus, series, mean_series);
    series = bsxfun(@times, series, 1./(std_series+eps) );
end

% Perform Kmeans
if(num_clusters > 1)
    % Normalize to zero mean across subjects
    series = bsxfun(@minus, series, mean(series, 2) );
    tic, clustered = direcClus_fix_bessel_bsxfun(series, num_clusters, size(series, 2) - 1, num_tries, 0, 0, 0, 1e-4, 1, max_iter, 1); toc
    cidx = clustered.clusters;
    lambda = clustered.lambda;
    mtc = clustered.mtc;
    p = clustered.p;
    lowerbound = clustered.likelihood(end);
    
    if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
        new_cidx = zeros(length(non_zero_corr_index), 1);
        new_cidx(non_zero_corr_index) = cidx;
        cidx = new_cidx; 
    end
else
    cidx = ones(size(series, 1), 1);
    lambda = 0;
    mtc = 0;
    p = 0;
end

% Write Clustering Results
lh_labels = zeros(lh_num_verts, 1);
lh_labels(l1) = cidx(1:length(l1));

rh_labels = zeros(rh_num_verts, 1);
rh_labels(l2) = cidx(length(l1)+1:end);

save(output_file, 'lh_labels', 'rh_labels', 'lambda', 'mtc', 'lowerbound');


if(num_clusters > 1)
    tic; s = silhouette(series, cidx(non_zero_corr_index), 'correlation'); toc
    if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
        new_s = ones(length(non_zero_corr_index), 1);
        new_s(non_zero_corr_index) = s;
        s = new_s; 
    end
else
    s = zeros(size(series, 1), 1);
end
lh_s = ones(lh_num_verts, 1);
lh_s(l1) = s(1:length(l1));

rh_s = ones(rh_num_verts, 1);
rh_s(l2) = s(length(l1)+1:end);

save(output_file, 'lh_labels', 'rh_labels', 'lh_s', 'rh_s', 'lambda', 'mtc', 'p', 'lowerbound');
