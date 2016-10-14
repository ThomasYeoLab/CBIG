function log_p = CBIG_IntersectALELogProb(ale_log_mat)

% log_p = CBIG_IntersectALELogProb(ale_log_mat)
%
% Compute a summation of log activation probabilities for a given set of brain locations
% FORMAT log_p = CBIG_IntersectALELogProb(ale_log_mat)
%
% ale_log_mat  = #foci    x # voxels/vertices 
%             or #studies x # voxels/vertices matrix of log activation probabilities
% 
% log_p        = 1        x # voxels/vertices vector of log activation probabilities
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  log_p = sum(ale_log_mat, 1); 
