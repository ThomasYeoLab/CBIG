function CBIG_VonmisesSeriesConsistencySurf(mesh_name, mask, num_clusters, output_dir, profile1, profile2, num_smooth, num_tries, rand_num, dim, normalize)

% CBIG_VonmisesSeriesConsistencySurf(mesh_name, mask, num_clusters, output_dir, profile1, profile2, num_smooth, num_tries, rand_num, dim, normalize)
% 
% Compute the consistency of von Mises-Fisher surface parcellation. The
% vertices in surface correlation profiles are splited into two branches.
% Clustering algorithm is performed on these two branches seperately. The
% consistency is indicated by the cost computed from Hungarian matching
% between the parcellation results from two branches.
%
% rand_num: number of differenct split of data
% dim:      1 or 2, whether differenct vertices are on the first dimesion
%           or the second dimension in input profiles
% 
% For further explanation of input arguments, please see
% CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun.m 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



a1 = strfind(profile1, 'lh.');
a2 = strfind(profile1, '.mgz');
out_str = profile1(a1+3:a2-1);

if(ischar(num_clusters))
    num_clusters = str2num(num_clusters);
end

if(num_clusters == 1)
    error('Does not handle single cluster');
end

if(ischar(normalize))
   normalize = str2num(normalize); 
end

if(ischar(num_tries))
    num_tries = str2num(num_tries);
end

if(ischar(num_smooth))
    num_smooth = str2num(num_smooth);
end

if(ischar(rand_num))
    rand_num = str2num(rand_num);
end

if(rand_num == 0)
    error(['Assumes rand num (' num2str(rand_num) ') greater than 0!']);
end

if(ischar(dim))
    dim = str2num(dim);
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
    
% znormalize (series assumed to be voxels x subjects or voxels x profile)
if(normalize)
    mean_series = nanmean(series, 1);
    std_series = nanstd(series, 1, 1);
    series = series - repmat(mean_series, size(series, 1), 1);
    series = series./repmat(std_series+eps, size(series, 1), 1);
end

% normal clustering
output_file = [output_dir '/Cluster' num2str(num_clusters, '%03d') '.s' num2str(num_smooth, '%02d') '.tries' num2str(num_tries) '.rand' num2str(rand_num, '%03d') '.znorm' num2str(normalize) '.dim' num2str(dim) '.' out_str '.mat'];
if(exist(output_file, 'file'))
    load(output_file);
    start = length(con_struct.orig_overlap) + 1;
else
    start = 1;
end

for i = start:rand_num

    rand('state', 100*i);
    disp(num2str(i)); % split data
    p = randperm(size(series, dim));
    idx1 = p(1:round(size(series, dim)/2));
    idx2 = p(round(size(series, dim)/2)+1:end);

    if(dim == 2)
        series1 = series(:, idx1);
        series2 = series(:, idx2);
    else
        series1 = series(idx1, :);
        series2 = series(idx2, :);
    end

    % cluster
    series1 = series1 - repmat(mean(series1, 2), 1, size(series1, 2));
    clustered1 = direcClus(series1, num_clusters, size(series1, 2) - 1, num_tries, 0, 0, 0, 1e-4, 1, 50, 0);

    series2 = series2 - repmat(mean(series2, 2), 1, size(series2, 2));
    clustered2 = direcClus(series2, num_clusters, size(series2, 2) - 1, num_tries, 0, 0, 0, 1e-4, 1, 50, 0);

    % output
    con_struct.out{i, 1}.idx = idx1;
    con_struct.out{i, 1}.clusters = clustered1.clusters;
    con_struct.out{i, 1}.p = clustered1.p;
    con_struct.out{i, 1}.mtc = clustered1.mtc;
    con_struct.out{i, 1}.lambda = clustered1.lambda;

    con_struct.out{i, 2}.idx = idx2;
    con_struct.out{i, 2}.clusters = clustered2.clusters;
    con_struct.out{i, 2}.p = clustered2.p;
    con_struct.out{i, 2}.mtc = clustered2.mtc;
    con_struct.out{i, 2}.lambda = clustered2.lambda;

    if(dim == 2)
        % match between two clusterings
        [labels, assign, cost] = CBIG_HungarianClusterMatch(clustered1.clusters, clustered2.clusters);
        con_struct.orig_overlap(i) = -cost/length(clustered1.clusters);
    else
        % Compute clustering on set 1 using parameters from set 2.
        series1 = series1 ./ repmat(sqrt(sum(series1.^2, 2)),1,size(series1, 2));
        distctrs1 = clustered2.lambda*series1*(clustered2.mtc');
        maxdists1 = max(distctrs1,[],2);

        rr1 = repmat(clustered2.p, size(distctrs1, 1), 1) .* exp(distctrs1 - maxdists1*ones(1, num_clusters));
        sumrr1 = sum(rr1,2);
        prob2 = rr1 ./ sumrr1(:,ones(num_clusters,1));

        [maxr2, assign2] = max(prob2,[],2);

        % match predicted and actual clustering
        [labels, assign, cost] = CBIG_HungarianClusterMatch(clustered1.clusters, assign2);
        con_struct.orig_overlap(i) = -cost/length(clustered1.clusters);
    end

    % compute random partition
    rand_clusters1 = randsample(num_clusters, length(clustered1.clusters), true);
    rand_clusters2 = randsample(num_clusters, length(clustered1.clusters), true);
    [labels, assign, cost] = CBIG_HungarianClusterMatch(rand_clusters1, rand_clusters2);
    con_struct.rand_overlap(i) = -cost/length(clustered1.clusters);

    save(output_file, 'con_struct');
end





