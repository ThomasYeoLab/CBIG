function G = CBIG_SPGrad_create_graph(sbjMesh, grad_data)

% G = CBIG_SPGrad_create_graph(sbjMesh, grad_data)
%
% This script creates a graph structure based on the given data and surface mesh
%
% Input:
%     - resolution: (scalar)
%       the number of vertices of surface mesh, e.g. 40962.
%
%     - output_dir:
%       the output directory which saves the cifti template.
%
% Output:
%     - G:
%       the constructed graph.
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

verts = 1:size(sbjMesh.vertexNbors,2);

for n = 1:size(sbjMesh.vertexNbors,1)
    neighbors = sbjMesh.vertexNbors(n,:);
    mask = neighbors ~= 0;
    grad_weight = (grad_data(mask) + grad_data(neighbors(mask)))./2;
    if(n==1)
        G = graph(verts(mask), neighbors(mask), grad_weight);
    else
        G = addedge(G, verts(mask), neighbors(mask), grad_weight);
    end
end
    

