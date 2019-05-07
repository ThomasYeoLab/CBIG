function CBIG_MMLDA_brain_to_zscore_create_refpara(vol4D_bl, mask, ind_CN_bl_brain, regressor_bl_brain, zscore_method, output)
% CBIG_MMLDA_brain_to_zscore_create_refpara(vol4D_bl, mask, ind_CN_bl_brain, regressor_bl_brain, zscore_method, output)
%
% Create a mat file which contains reference parameters for CBIG_MMLDA_brain_to_doc.m function.
% The variable [ref_mean] [ref_beta] will be used for regression with respect to reference cohort.
% The variable [ref_reg_mean] [ref_reg_std] will be used for z score with respect to reference cohort. 
%
% Input:
%   - vol4D_bl              : mri structure (from MRIread) of MRI 4D baseline volume 
%   - mask                  : mri structure (from MRIread) of Grey Matter mask
%   - ind_CN_bl_brain       : index of cognitive normal group at baseline brain modality and we use CN as reference cohort. 
%   - regressor_bl_brain    : N x M matrix, N is # of subjects at baseline brain modality
%   - zscore_method         : 'meanCNstdCN' or 'meanCNstdALL'. Use mean of cognitive normal (CN) and stand deviation of 
%                             cognitive normal (CN) or all subjects (ALL) to do the z score.
%   - output                : full path of a .mat file (e.g. ~/example/brain_ref.mat)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Load bl 4D MRI
vols_bl = vol4D_bl.vol;
vol_size = size(vols_bl);
% Load 3D mask
mask_3d = mask.vol;
mask_size = size(mask_3d);
mask_3d(mask_3d~=0) = 1; %binarize the mask, spm mask have values 0.999 
if vol_size(1:3)~=mask_size
    error('Error: mask should have the same size as input volume.')
end

% Mask the brain vol
vol_2d = reshape(vols_bl,[size(vols_bl,1)*size(vols_bl,2)*size(vols_bl,3) size(vols_bl,4)]);
vol_2d = vol_2d';
mask_1d = mask_3d(:);
mask_1d = mask_1d';
sub_vols_bl = vol_2d(:,logical(mask_1d));
sub_vols_bl(sub_vols_bl==0) = 1; % avoid log10() to be infinity

% apply log10 
sub_vols_bl = log(sub_vols_bl)/log(10);

%%%
% Regress the nuisance regressor
%%%
% Estimate GLM parameters with only CN subjects
Y_CN = sub_vols_bl(ind_CN_bl_brain,:);
ref_mean = mean(Y_CN,1);
X_CN = [ones(size(Y_CN,1),1) regressor_bl_brain(ind_CN_bl_brain,:)];
ref_beta = (X_CN'*X_CN)\(X_CN'*Y_CN);

% Regress out regressors for all subjects by computing Y - X*ref_beta +
% ref_mean for each voxel
Y = sub_vols_bl;
X = [ones(size(Y,1),1) regressor_bl_brain];
sub_vols_reg_bl = bsxfun(@plus, Y-X*ref_beta, ref_mean);

%%%
% z score with respect to CN
%%%
if strcmp(zscore_method, 'meanCNstdCN') 
    % compute post-regression CN group mean and std for z-normalization
    ref_reg_mean = mean(sub_vols_reg_bl(ind_CN_bl_brain,:),1);
    ref_reg_std = std(sub_vols_reg_bl(ind_CN_bl_brain,:),0,1);
elseif strcmp(zscore_method, 'meanCNstdALL')
    % minus mean(CN) divide std(ALL)
    ref_reg_mean = mean(sub_vols_reg_bl(ind_CN_bl_brain,:),1);
    ref_reg_std = std(sub_vols_reg_bl,0,1);
else
    error('Error: No such zscore_method option.')
end

%%%
% save reference parameters into a mat file
%%%
save(output, 'ref_mean', 'ref_beta', 'ref_reg_mean', 'ref_reg_std')