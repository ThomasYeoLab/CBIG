function [corr_lh,corr_rh,corr_lh_ex,corr_rh_ex] = CBIG_preproc_compare_two_surfaces(lh_surface_file1,rh_surface_file1,lh_surface_file2,rh_surface_file2,fs_dim)

% [corr_lh,corr_rh,corr_lh_ex,corr_rh_ex] = CBIG_preproc_compare_two_surfaces(lh_surface_file1,rh_surface_file1,lh_surface_file2,rh_surface_file2)
%
% Given two surface file name (lh & rh), calculate the Pearson correlation between each corresponding voxel.
% Input:
%   - lh_surface_file1: file name of 1st surface left hemisphere  
%   - rh_surface_file1: file name of 1st surface right hemisphere  
%   - lh_surface_file2: file name of 2nd surface left hemisphere  
%   - rh_surface_file2: file name of 2nd surface right hemisphere 
%   - fs_dim: fsaverage4, fsaverage5 or fsaverage6
% Output:
%   - corr_lh: correlation of 2 surface lh (set the corr of medial wall as 0); a 1xN matrix, N is num of voxels in left hemisphere
%   - corr_rh: correlation of 2 surface rh (set the corr of medial wall as 0); a 1xN matrix, N is num of voxels in right hemisphere
%   - corr_lh_ex: correlation of 2 surface lh (exclude medial wall); a 1xM matrix, M is num of voxels in left hemisphere excluding the medial wall
%   - corr_rh_ex: correlation of 2 surface rh (exclude medial wall); a 1xM matrix, M is num of voxels in right hemisphere excluding the medial wall
%
% Author: Nanbo Sun
% Date: 2016/06/15


%% load surface1
lh_fmri = MRIread(lh_surface_file1);
rh_fmri = MRIread(rh_surface_file1);

lh_vol = lh_fmri.vol;
rh_vol = rh_fmri.vol;

vol_size = size(lh_vol);

lh_vol2d = reshape(lh_vol,[vol_size(1)*vol_size(2)*vol_size(3) vol_size(4)]);
rh_vol2d = reshape(rh_vol,[vol_size(1)*vol_size(2)*vol_size(3) vol_size(4)]);

lh_vol2d_1 = lh_vol2d';
rh_vol2d_1 = rh_vol2d';


%% load surface2
lh_fmri = MRIread(lh_surface_file2);
rh_fmri = MRIread(rh_surface_file2);

lh_vol = lh_fmri.vol;
rh_vol = rh_fmri.vol;

vol_size = size(lh_vol);

lh_vol2d = reshape(lh_vol,[vol_size(1)*vol_size(2)*vol_size(3) vol_size(4)]);
rh_vol2d = reshape(rh_vol,[vol_size(1)*vol_size(2)*vol_size(3) vol_size(4)]);

lh_vol2d_2 = lh_vol2d';
rh_vol2d_2 = rh_vol2d';

%% calculate correlation
% read medial wall label
lh_mesh = CBIG_ReadNCAvgMesh('lh',fs_dim,'inflated','cortex');
rh_mesh = CBIG_ReadNCAvgMesh('rh',fs_dim,'inflated','cortex');

lh_medial_wall_label = lh_mesh.MARS_label;
rh_medial_wall_label = rh_mesh.MARS_label;

% exclude the medial wall
lh_vol2d_1_ex = lh_vol2d_1(:,lh_medial_wall_label==2);
rh_vol2d_1_ex = rh_vol2d_1(:,rh_medial_wall_label==2);

lh_vol2d_2_ex = lh_vol2d_2(:,lh_medial_wall_label==2);
rh_vol2d_2_ex = rh_vol2d_2(:,rh_medial_wall_label==2);

corr_lh_ex = CBIG_preproc_corr_matrix(lh_vol2d_1_ex,lh_vol2d_2_ex);
corr_rh_ex = CBIG_preproc_corr_matrix(rh_vol2d_1_ex,rh_vol2d_2_ex);

if sum(isnan(corr_lh_ex)) > 0
    fprintf('Warning: lh correlation have NaNs\n')
elseif sum(isnan(corr_rh_ex)) > 0
    fprintf('Warning: rh correlation have NaNs\n')
end

% set the medial wall corr as 0
corr_lh = zeros(size(lh_medial_wall_label));
corr_rh = zeros(size(rh_medial_wall_label));

corr_lh(lh_medial_wall_label==2) = corr_lh_ex;
corr_rh(rh_medial_wall_label==2) = corr_rh_ex;

corr_lh_ex = corr_lh_ex(~isnan(corr_lh_ex));
corr_rh_ex = corr_rh_ex(~isnan(corr_rh_ex));

% corr_lh=corr_matrix(lh_vol2d_1,lh_vol2d_2);
% corr_rh=corr_matrix(rh_vol2d_1,rh_vol2d_2);
% % exclude the NaN value
% corr_lh_ex=corr_lh(~isnan(corr_lh));
% corr_rh_ex=corr_rh(~isnan(corr_rh));
% 
% % set the NaN value as 0
% corr_lh(isnan(corr_lh))=0;
% corr_rh(isnan(corr_rh))=0;











