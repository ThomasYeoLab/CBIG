function null_p = CBIG_GenerateALENullDist_EickhoffUnion(ale_mat, num_samples)

% null_p = CBIG_GenerateALENullDist_EickhoffUnion(ale_mat, num_samples)
%
% Generate a union of activation probabilities that were sampled from 
% a null hypothesis of spatial independence
% FORMAT null_p = CBIG_GenerateALENullDistr_EickhoffUnion(ale_mat, num_samples)
%
% ale_mat = # studies x # relevant voxels matrix of activation probabilities
% num_samples = # sampled locations
% 
% null_p = 1 x num_samples vector of union of respective acitvation probabilities
%          for the sampled locations
% ____________________________________________________________________________________________________
%
% CBIG_GenerateALENullDist_EickhoffUnion returns a union of activation probabilities, which were sampled from
% random, spatially independent brain locations as described in Eickhoff et. al, 2009.
% ____________________________________________________________________________________________________
%
% Refs:
%
% Eickhoff et al., Human Brain Mapping, 2009. Coordinate-based activation likelihood estimation meta-analysis
% of neuroimaging data: A random-effects approach based on empirical estimates of spatial uncertainty
%
% Written by Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_studies = size(ale_mat, 1);
num_voxels  = size(ale_mat, 2);

samples     = randi(num_voxels, [num_studies num_samples]);
study_index = repmat(transpose(1:num_studies), [1 num_samples]);
samples     = sub2ind(size(ale_mat), study_index(:), samples(:));

null_p = reshape(ale_mat(samples), [num_studies num_samples]);
null_p = 1 - prod(1 - null_p, 1);
