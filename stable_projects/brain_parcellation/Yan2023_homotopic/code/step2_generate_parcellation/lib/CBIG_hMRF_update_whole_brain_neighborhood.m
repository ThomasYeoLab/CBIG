function Neighborhood = CBIG_hMRF_update_whole_brain_neighborhood(Neighborhood_lh, Neighborhood_rh,...
     Neighborhood_lh_rh, c, d)
% Neighborhood = CBIG_hMRF_update_whole_brain_neighborhood(Neighborhood_lh, Neighborhood_rh,...
% Neighborhood_lh_rh, c, d)
%
% This function constructs a sparse matrix with parameters c and d based on weighed neighborhood matrices.

% For the notations below:
% NL = no of vertices per hemisphere;
% NR = no of vertices per hemisphere;

% Input
%   - Neighborhood_lh: (matrix)
%     NL x NL matrix representing the left hemisphere neighborhood.
%   - Neighborhood_rh: (matrix)
%     NR x NR matrix representing the right hemisphere neighborhood.
%   - Neighborhood_lh_rh: (matrix)
%     NL x NR matrix representing the interhemisphere connections.
%   - c: (double)
%     The hyperparameter for the gradient-weighted MRF term in the cost function.
%     The larger the value, the stronger the gradient-weighted MRF term would weigh in the cost function, thus
%     the parcellation tend to conform more to the gradient boundaries.
%   - d: (double)
%     The hyperparameter that controls interhemispheric symmetry in the cost function.
%     The stronger d is, the stronger the spatial symmetry constraint is (i.e., resultant parcellation more symmetric).
% Output
%   - Neighborhood: (matrix)
%     (NL+NR) x (NL+NR) matrix representing the neighborhood for both hemispheres,
%      both interhemispheric and intrahemispheric.
%
% Example
%   - Neighborhood = CBIG_hMRF_update_whole_brain_neighborhood(Neighborhood_lh,...
%     Neighborhood_rh, Neighborhood_lh_rh, c, d)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

Neighborhood = blkdiag(Neighborhood_lh, Neighborhood_rh) .* (c/d);

lh_cortex_vertices = size(Neighborhood_lh_rh, 1);
Neighborhood(1:lh_cortex_vertices, lh_cortex_vertices+1:end) = Neighborhood_lh_rh; 
Neighborhood(lh_cortex_vertices+1:end, 1: lh_cortex_vertices) = Neighborhood_lh_rh';

Neighborhood = triu(Neighborhood);
Neighborhood = sparse(Neighborhood);
end