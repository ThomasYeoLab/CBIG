function null_log_p = CBIG_GenerateALENullDist_EickhoffIntersect(ale_log_mat, num_samples)

% null_log_p = CBIG_GenerateALENullDist_EickhoffIntersect(ale_log_mat, num_samples)
%
% Generate a sum of log activation probabilities that were sampled from 
% a null hypothesis of spatial independence
% FORMAT null_p = CBIG_GenerateALENullDistr_EickhoffIntersect(ale_log_mat, num_samples)
%
% ale_log_mat = # studies x # relevant voxels matrix of log activation probabilities
% num_samples = # sampled locations
% 
% null_p = 1 x num_samples vector of sum of respective log acitvation probabilities
%          for the sampled locations
% ____________________________________________________________________________________________________
%
% CBIG_GenerateALENullDist_EickhoffIntersect returns a sum of log activation probabilities, which were sampled from
% random, spatially independent brain locations as described in Eickhoff et. al, 2009.
% ____________________________________________________________________________________________________
%
% Refs:
%
% Eickhoff et al., Human Brain Mapping, 2009. Coordinate-based activation likelihood estimation meta-analysis
% of neuroimaging data: A random-effects approach based on empirical estimates of spatial uncertainty
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_studies = size(ale_log_mat, 1);
num_voxels  = size(ale_log_mat, 2);

samples     = randi(num_voxels, [num_studies num_samples]);
study_index = repmat(transpose(1:num_studies), [1 num_samples]);
samples     = sub2ind(size(ale_log_mat), study_index(:), samples(:));

null_log_p = reshape(ale_log_mat(samples), [num_studies num_samples]);
null_log_p = sum(null_log_p, 1);
