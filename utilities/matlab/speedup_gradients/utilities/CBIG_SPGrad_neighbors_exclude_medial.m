function reorg_neighbors = CBIG_SPGrad_neighbors_exclude_medial(orig_neighbors, mask)

% reorg_neighbors = CBIG_SPGrad_neighbors_exclude_medial(orig_neighbors, mask)
%
% This script excludes medial wall vertices from the neighborhood structure and make
% the neighborhood vertices indices contiguous.
%
% Input:
%     - orig_neighbors:
%       the original neighborhood structure.
%
%     - mask: (#num_vertices x 1 binary vector)
%       the medial area mask. A #num_vertices x 1 binary vector. The medial area 
%       should be denoted as 1, others should be denoted as 0.
%
% Output:
%     - reorg_neighbors:
%       the neighborhood structure excluding medial wall vertices
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%orig_neighbors: KxN
%mask: Nx1

medial_set = find(mask==1);
orig_neighbors(:,mask) = [];
%medial_neigh: KxN logical variable, 1 indicates medial wall vertices
medial_neigh = ismember(orig_neighbors, medial_set);
%if medial wall vertices are neighbors of other vertices
orig_neighbors(medial_neigh) =  NaN;


%make indices contiguous
reorg_neighbors = orig_neighbors;
for i=1:length(medial_set)
    reorg_neighbors(orig_neighbors > medial_set(i)) = reorg_neighbors(orig_neighbors > medial_set(i)) - 1;
end


end