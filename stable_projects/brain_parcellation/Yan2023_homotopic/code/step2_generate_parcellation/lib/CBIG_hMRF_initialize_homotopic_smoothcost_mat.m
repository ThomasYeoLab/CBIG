function Smooth_cost = CBIG_hMRF_initialize_homotopic_smoothcost_mat(num_cluster_per_hemi, d)
% Smooth_cost = CBIG_hMRF_initialize_homotopic_smoothcost_mat(num_cluster_per_hemi, d)
%
% This function initializes smoothcost matrix with 0 entries for the main diagonal 
% and sub-diagonals between hemispheres.
% For example, if lh and rh have 3 clusters respectively, the smoothcost matrix is given by:
% 0 d d 0 d d
% d 0 d d 0 d
% d d 0 d d 0
% 0 d d 0 d d
% d 0 d d 0 d
% d d 0 d d 0
% 

% For the notations below:
% k = total no of parcels 

% Input
%   - num_cluster_per_hemi: (integer)
%     The number of clusters per hemisphere.
%   - d: (double)
%     The hyperparameter d controlling homotopy.
% Output
%   - Smooth_cost: (matrix)
%     A kxk matrix representing smoothcost and to be fed to the GCO algorithm.
%
% Example
%   - Smooth_cost = CBIG_hMRF_initialize_homotopic_smoothcost_mat(200, 1000)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_cluster = 2 * num_cluster_per_hemi;
Smooth_cost = ones(num_cluster, num_cluster, 'single') * single(d); % sets the diagonal to 0
Smooth_cost(1: (num_cluster+1): end) = 0; % set the diagonals of the cross-hemisphere matrices to zero
Smooth_cost(num_cluster_per_hemi+1: num_cluster+1: num_cluster * num_cluster_per_hemi) = 0; 
Smooth_cost(num_cluster * num_cluster_per_hemi+1: num_cluster+1: end) = 0;
end
