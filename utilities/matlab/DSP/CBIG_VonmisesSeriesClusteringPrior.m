function CBIG_VonmisesSeriesClusteringPrior(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, normalize, prior_type, input_surf_clusters, max_iter)

% CBIG_VonmisesSeriesClusteringPrior(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, normalize, prior_type, input_surf_clusters, max_iter)
% 
% Von Mises-Fisher clustering with prior information.
% 
% prior_type:          "full" or "init"
%                      If priot_type = "full", this function inherits prior clustering results;
%                      If prior_type = "init", the prior information is used as starting point for further clustering.
% input_surf_clusters: surface clustering prior
%
% For further explanation of input arguments, please see CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun.m
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% prior_type = 'full' or 'init'

if(~exist('max_iter', 'var'))
   max_iter = 1000; % note increase in max-iteration
else
   if(ischar(max_iter))
       max_iter = str2num(max_iter);
   end
end

if(~strcmp(prior_type, 'full') && ~strcmp(prior_type, 'init'))
   error('Prior_type needs to be full or init'); 
end

% We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
rand('twister', 5489)

if(ischar(normalize))
   normalize = str2num(normalize); 
end

if(ischar(num_clusters))
   num_clusters = str2num(num_clusters); 
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
    series = series - repmat(mean_series, size(series, 1), 1);
    series = series./repmat(std_series+eps, size(series, 1), 1);
end

if(num_clusters > 1)
    % read surface clusters
    surf_clusters = load(input_surf_clusters);
    
    if(strcmp(prior_type, 'full'))
        disp('full');
        
        % von_mises pre-filtering
        series = series - repmat(mean(series, 2), 1, size(series, 2));
        series = series ./ repmat(sqrt(sum(series.^2, 2)),1,size(series, 2));
        
        % classify brain volume
        num_clusters = size(surf_clusters.mtc, 1);
        distctrs = surf_clusters.lambda*series*(surf_clusters.mtc');
        maxdists = max(distctrs,[],2);
        
        rr = exp(distctrs - maxdists*ones(1, num_clusters));
        sumrr = sum(rr,2);
        prob = rr ./ sumrr(:,ones(num_clusters,1));
        
        [tmp, cidx] = max(prob,[],2);
        lambda = surf_clusters.lambda; %inherit
        mtc = surf_clusters.mtc;
        p = surf_clusters.p;
    elseif(strcmp(prior_type, 'init'))
 
        disp('init');
        num_tries = 1;
        series = series - repmat(mean(series, 2), 1, size(series, 2));
        tic, clustered = direcClus_fix_bessel_bsxfun(series, num_clusters, size(series, 2) - 1, num_tries, surf_clusters.lambda, surf_clusters.p, surf_clusters.mtc, 1e-6, 1, max_iter, 1); toc
        cidx = clustered.clusters;
        lambda = clustered.lambda;
        mtc = clustered.mtc;
        p = clustered.p;
    end
else
    cidx = ones(size(series, 1), 1);
    lambda = 0;
    mtc = 0;
    p = 0;
end

if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
    new_cidx = zeros(length(non_zero_corr_index), 1);
    new_cidx(non_zero_corr_index) = cidx;
    cidx = new_cidx;
end    


% Write Clustering Results
lh_labels = zeros(lh_num_verts, 1);
lh_labels(l1) = cidx(1:length(l1));

rh_labels = zeros(rh_num_verts, 1);
rh_labels(l2) = cidx(length(l1)+1:end);

save(output_file, 'lh_labels', 'rh_labels', 'lambda', 'mtc', 'p');

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

save(output_file, 'lh_labels', 'rh_labels', 'lh_s', 'rh_s', 'lambda', 'mtc', 'p');
