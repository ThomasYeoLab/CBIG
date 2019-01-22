function CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, num_tries, normalize, max_iter, no_silhouette)

% CBIG_VonmisesSeriesClustering_fix_bessel_randnum_bsxfun(mesh_name, mask, num_clusters, output_file, profile1, profile2, num_smooth, num_tries, normalize, max_iter)
%
% Von Mises-Fisher clustering on surface data.
% 
% Input arguments:
%     - mesh_name     : e.g. 'fsaverage5'
%     - mask          : e.g. 'cortex' if mesh_name is 'fsaverage*';
%                       otherwise, pass in empty string or 'NONE'
%     - num_clusters  : number of clusters
%     - output_file   : output file name
%     - profile1      : group average profile on left hemisphere for data in
%                       fsaverage* space; 
%                       or group average profile on entire
%                       cortex for data in fs_LR_32k space.
%     - profile2      : group average profile on right hemisphere for data
%                       in fsaverage* space;
%                       it is not useful for data in fs_LR_32k space. You
%                       can pass in empty string or 'NONE'.
%     - num_smooth    : how many times that smoothing is performed. It is
%                       only used when mesh_name is 'fsaverage*'
%     - num_tries     : number of difference random initialization
%     - normailize    : 0 or 1, whether z-normalization is performed across vertices
%     - max_iter      : maximum number of iterations for one random initialization
%     - no_silhouette : if 1 is passed in for this flag, silhouette step
%                       will be skipped. Default is 0, meaning silhouette
%                       step will be performed when num_clusters>1 and data
%                       are in fsaverage* space.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~isempty(strfind(mesh_name, 'fsaverage')))
    lambda = 500;            % function direcClus_fix_bessel_bxfun() treats 0 as not specified, lambda will be set as its default (500)
elseif(strcmp(mesh_name, 'fs_LR_32k'))
    lambda = 650;          % For data in fs_LR_32k space, lambda is set to be 650
else
    error('Unknown mesh name.')
end


if(~exist('max_iter', 'var'))
   max_iter = 100; 
else
   if(ischar(max_iter))
       max_iter = str2num(max_iter);
   end
end

if(~exist('no_silhouette', 'var'))
    no_silhouette = 0;     % silouette step will be performed 
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
if(~isempty(strfind(mesh_name, 'fsaverage')))
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_name, 'inflated', mask);
    l1 = find(lh_avg_mesh.MARS_label == 2); lh_num_verts = size(lh_avg_mesh.vertices, 2);
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_name, 'inflated', mask);
    l2 = find(rh_avg_mesh.MARS_label == 2); rh_num_verts = size(rh_avg_mesh.vertices, 2);
    l = [l1 l2+length(lh_avg_mesh.MARS_label)];
else
    % mesh.l: 0 - medial wall, 1 - cortex
    [mesh.v, mesh.l, mesh.ct] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'medialwall.annot'));
    lh_num_verts = length(mesh.l) / 2;
    rh_num_verts = lh_num_verts;
    cort_label = mesh.ct.table(2, 5);
    l1 = 1:lh_num_verts;          l1 = l1(mesh.l(l1)==cort_label);
    l2 = 1:rh_num_verts;          l2 = l2(mesh.l(l2+lh_num_verts)==cort_label);
    l = [l1 l2+lh_num_verts];
end

% read data (voxels x N subjects)
series = read_fmri(profile1);
if(~isempty(strfind(mesh_name, 'fsaverage')))
    series2 = read_fmri(profile2);
    series = [series; series2];
end


% We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
% rand('twister',5489) was moved to after MRIread is because MRIread calls/apps/arch/Linux_x86_64/freesurfer/4.5.0/matlab/load_nifti.m, which calls rand('state', sum(100*clock));
rand('twister',5489)

% smooth, only applied for data in fsaverage* space
if(~isempty(strfind(mesh_name, 'fsaverage')))
    series(1:end/2, :)     = transpose(MARS_AverageData(lh_avg_mesh, transpose(series(1:end/2, :)), 0, num_smooth));
    series(end/2+1:end, :) = transpose(MARS_AverageData(rh_avg_mesh, transpose(series(end/2+1:end, :)), 0, num_smooth));
end

% extract mask voxels series
series = series(l, :);

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
    tic, clustered = direcClus_fix_bessel_bsxfun(series, num_clusters, size(series, 2) - 1, num_tries, lambda, 0, 0, 1e-4, 1, max_iter, 1); toc
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
    cidx = ones(length(l), 1);
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


if(num_clusters > 1 && ~isempty(strfind(mesh_name, 'fsaverage')) && no_silhouette==0)
    tic; s = silhouette(series, cidx(non_zero_corr_index), 'correlation'); toc
    if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
        new_s = ones(length(non_zero_corr_index), 1);
        new_s(non_zero_corr_index) = s;
        s = new_s; 
    end
else
    s = zeros(length(l), 1);
end
lh_s = ones(lh_num_verts, 1);
lh_s(l1) = s(1:length(l1));

rh_s = ones(rh_num_verts, 1);
rh_s(l2) = s(length(l1)+1:end);

save(output_file, 'lh_labels', 'rh_labels', 'lh_s', 'rh_s', 'lambda', 'mtc', 'p', 'lowerbound');




function vol = read_fmri(fmri_name)

% [fmri, vol] = read_fmri(fmri_name)
% Given the name of averaged correlation profile file (fmri_name), this
% function read in the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%       It can be .nii.gz file for correlation profile in fsaverage*
%       spaces; or .mat file for correlation profile in fs_LR_32k space
%       (there must be a variable called 'profile_mat' in the .mat file).
%
% Output:
%     - vol:
%       A num_voxels x num_timepoints matrix which is the reshaped 'vol'
%       structure for .nii.gz file or the variable 'profile_mat' for .mat file.
%

if (~isempty(strfind(fmri_name, '.nii.gz')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = fmri.vol;
    vol_size = size(vol);
    if(length(vol_size) < 4)
        vol = reshape(vol, prod(vol_size(1:3)), 1);
    else
        vol = reshape(vol, prod(vol_size(1:3)), vol_size(4));
    end
    fmri.vol = [];
else
    load(fmri_name);
    vol = profile_mat;
end
