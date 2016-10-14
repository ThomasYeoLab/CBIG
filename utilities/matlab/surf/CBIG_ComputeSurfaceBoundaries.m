function BoundaryVec = CBIG_ComputeSurfaceBoundaries(mesh, labels)

% BoundaryVec = CBIG_ComputeSurfaceBoundaries(mesh, labels)
% This function is used to generate a two-vertex thick boundary defined by
% labels.
% Example:
% lh_avg_mesh= CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
% BoundaryVec = CBIG_ComputeSurfaceBoundaries(lh_avg_mesh,lh_labels);
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


BoundaryVec = zeros(length(labels), 1);
maxNeighbors = size(mesh.vertexNbors, 1);

for i = 1:length(labels)
    
    label_vertex = int32(labels(i));
    
    for k = 1:maxNeighbors
        v_neighbor = mesh.vertexNbors(k, i);
        if(v_neighbor ~= 0 && int32(labels(v_neighbor)) ~= label_vertex)
            BoundaryVec(i) = 1;
        end
    end
end

BoundaryVec = logical(BoundaryVec);