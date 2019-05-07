function Z_flip = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara)
% Z_flip = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara)
%
% Convert brain 4D volume to zscore with respect to reference cohort.
% 1. Take log10
% 1. Regress out regressors with respect to reference.
% 2. Zscore withe respect ot reference.
% 3. Flip z score so that larger means worse atrophy 
%
% Input:
%   - vol4D                 : mri structure (from MRIread) of MRI 4D volume, which is a concatenation of 3D volumes.
%   - mask                  : mri structure (from MRIread) of Grey Matter mask
%   - regressor             : N x M matrix, N is # of subjects, M is # of regressors
%   - refpara               : reference parameters from CBIG_MMLDA_brain_to_zscore_create_refpara.m
%
% Output:
%   - Z_flip                : N x A marix, N is # of subjects, A is # of voxels within the mask.
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% load viscode 4D MRI
vols = vol4D.vol;
vol_size = size(vols);
% Load 3D mask
mask_3d = mask.vol;
mask_size = size(mask_3d);
mask_3d(mask_3d~=0) = 1; %binarize the mask, spm mask have values 0.999 
if vol_size(1:3)~=mask_size
    error('Error: mask should have the same size as input volume.')
end

% Mask the brain vol
vol_2d = reshape(vols,[vol_size(1)*vol_size(2)*vol_size(3) vol_size(4)]);
vol_2d = vol_2d';
mask_1d = mask_3d(:);
mask_1d = mask_1d';
sub_vols = vol_2d(:,logical(mask_1d));
sub_vols(sub_vols==0) = 1; % avoid log10() to be infinity

% apply log10 
sub_vols = log(sub_vols)/log(10);

% load reference parameters
load(refpara)
%%%
% Regress the nuisance regressor
%%%
% Regress out regressors for all subjects by computing Y - X*b_CN +
% mean(Y_CN) for each voxel
Y = sub_vols;
X = [ones(size(Y,1),1) regressor];
sub_vols_reg = bsxfun(@plus, Y-X*ref_beta, ref_mean);

%%%
% z score with respect to reference 
%%%
Z = bsxfun(@minus, sub_vols_reg, ref_reg_mean);
Z = bsxfun(@rdivide, Z, ref_reg_std);

%%%
% Flip z score so that larger z score corresponds to more atrophy
%%%
Z_flip = -Z;

% Confirm that there is no NaN, Inf, -Inf in your doc
if find(isnan(Z))
    error('Error: Find NaN in Z score.\n')
end
if find(isinf(Z))
    error('Error: Find Inf or -Inf in Z score.\n')
end
