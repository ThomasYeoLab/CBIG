function [corr_vol,vol_size] = CBIG_preproc_compare_two_vols(vol_file1,vol_file2)

% [corr_vol,vol_size] = CBIG_preproc_compare_two_vols(vol_file1,vol_file2)
%
% Given two vol file name, calculate the Pearson correlation between each
% corresponding time course.
% Input:
%   - vol_file1: file name of 1st volume 
%   - vol_file2: file name of 2st volume  
% Output:
%   - corr_vol: Pearson correlation between each corresponding time course;
%   a 1xN matrix, N is num of voxels in vol
%   - vol_size: the size of the corr_vol
%
% Author: Nanbo Sun
% Date: 2016/06/15

%% load vol1
mri=MRIread(vol_file1);
mri_vol=mri.vol;
clear mri
mri_size=size(mri_vol);
vol1=reshape(mri_vol,[mri_size(1)*mri_size(2)*mri_size(3) mri_size(4)]);
clear mri_vol
vol1=single(vol1);
vol1=vol1';
fprintf('Load vol1 MNI done.\n')

%% load vol2
mri=MRIread(vol_file2);
mri_vol=mri.vol;
clear mri
mri_size=size(mri_vol);
vol2=reshape(mri_vol,[mri_size(1)*mri_size(2)*mri_size(3) mri_size(4)]);
clear mri_vol
vol2=single(vol2);
vol2=vol2';
fprintf('Load vol2 MNI done.\n')

%% calculate the correlation between two vols
corr_vol = CBIG_preproc_corr_matrix(vol1,vol2);
fprintf('corr done.\n')
vol_size = mri_size;



