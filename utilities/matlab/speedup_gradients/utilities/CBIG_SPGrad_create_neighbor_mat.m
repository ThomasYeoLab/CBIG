function neigbor_mat = CBIG_SPGrad_create_neighbor_mat(sbjMesh, grad_data)

% neigbor_mat = CBIG_SPGrad_create_neighbor_mat(sbjMesh, grad_data)
%
% This script creates a gradient neighborhood matrix for a given surface mesh
%
% Input:
%     - sbjMesh:
%       the surface mesh structure.
%
%     - grad_data:
%       the gradient map.
%
% Output:
%     - neighbor_mat:
%       the gradient neighborhood matrix. The egde weight between two nearby vertices is
%       the mean gradient value between these two vertices.
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
verts = 1:size(sbjMesh.vertexNbors,2);

for n = 1:size(sbjMesh.vertexNbors,1)
    neighbors = sbjMesh.vertexNbors(n,:);
    mask = neighbors ~= 0;
    grad_weight = (grad_data(mask) + grad_data(neighbors(mask)))./2;
    neighbor_mat(n,mask) = grad_weight;
end



    