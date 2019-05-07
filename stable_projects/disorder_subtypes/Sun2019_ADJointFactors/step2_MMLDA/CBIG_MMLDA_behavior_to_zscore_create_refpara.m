function CBIG_MMLDA_behavior_to_zscore_create_refpara(behavior_bl, ind_CN_bl_behavior, regressor_bl_behavior, zscore_method, output)
% CBIG_MMLDA_behavior_to_zscore_create_refpara(behavior_bl, ind_CN_bl_behavior, regressor_bl_behavior, zscore_method, output)
%
% Create a mat file which contains reference parameters for CBIG_MMLDA_behavior_to_doc.m function.
% The variable [ref_mean] [ref_beta] will be used for regression with respect to reference cohort.
% The variable [ref_reg_mean] [ref_reg_std] will be used for z score with respect to reference cohort. 
%
% Input:
%   - behavior_bl           : N x B matrix, N is # of subjects. B is # of behaviors (e.g. 28 in our paper). 
%                             We use all baseline subjects to compute the stand deviation, but you can use other subjects.
%   - ind_CN_bl_behavior    : N x 1 logical vector. Index of cognitive normal group at baseline behavior modality and 
%                             we use CN as reference cohort. 
%   - regressor_bl_behavior : N x M matrix, N is # of subjects at baseline behavior modality, M is # of regressors
%   - zscore_method         : 'meanCNstdCN' or 'meanCNstdALL'. Use mean of cognitive normal (CN) and stand deviation of 
%                             cognitive normal (CN) or all subjects (ALL) to do the z score.
%   - output                : full path of a .mat file 
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%%%
% regress out nuisance regressors
%%%
% Estimate GLM parameters with only CN subjects
Y_bl_CN_score = behavior_bl(ind_CN_bl_behavior, :);
ref_mean = mean(Y_bl_CN_score);
X_bl_CN = [ones(size(Y_bl_CN_score, 1), 1) regressor_bl_behavior(ind_CN_bl_behavior, :)];
ref_beta = (X_bl_CN'*X_bl_CN)\(X_bl_CN'*Y_bl_CN_score);

% Regress out regressors for all subjects by computing Y - X*b_CN +
% mean(Y_CN) for each voxel
Y = behavior_bl;
X = [ones(size(Y, 1), 1) regressor_bl_behavior];
behavior_reg_bl = bsxfun(@plus, Y-X*ref_beta, ref_mean);

%%%
% z score with respect to CN
%%%
if strcmp(zscore_method, 'meanCNstdCN')
    ref_reg_mean = mean(behavior_reg_bl(ind_CN_bl_behavior, :), 1);
    ref_reg_std = std(behavior_reg_bl(ind_CN_bl_behavior, :), 0, 1);
elseif strcmp(zscore_method, 'meanCNstdALL')
    ref_reg_mean = mean(behavior_reg_bl(ind_CN_bl_behavior, :), 1);
    ref_reg_std = std(behavior_reg_bl, 0, 1);
else
    error('Error: No such zscore_method option.')
end

%%%
% save reference parameters into a mat file
%%%
save(output, 'ref_mean', 'ref_beta', 'ref_reg_mean', 'ref_reg_std')



