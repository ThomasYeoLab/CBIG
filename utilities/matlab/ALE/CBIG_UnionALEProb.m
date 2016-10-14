function p = CBIG_UnionALEProb(ale_mat)

% p = CBIG_UnionALEProb(ale_mat)
%
% Compute a union of activation probabilities for a given set of brain locations
% FORMAT p = CBIG_UnionALEProb(ale_mat)
%
% ale_mat  = #foci    x # voxels/vertices 
%         or #studies x # voxels/vertices matrix of activation probabilities
% 
% p        = 1        x # voxels/vertices vector of activation probabilities
% 
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  p = 1 - prod(1 - ale_mat, 1);  
