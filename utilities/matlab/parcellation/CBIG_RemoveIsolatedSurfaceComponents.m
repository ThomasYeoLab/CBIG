function labels = CBIG_RemoveIsolatedSurfaceComponents(avg_mesh, labels, abs_threshold)

% Remove isolated surface components.
% 
%   labels = CBIG_RemoveIsolatedSurfaceComponents(avg_mesh, labels, abs_threshold)
%   Input:
%       avg_mesh        : average mesh read from CBIG_ReadNCAvgMesh.m
%       labels          : parcellation results
%       abs_threshold   : threshold (default 5)
%   Output:
%       labels          : return labels after removing
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(nargin < 3)
   abs_threshold = 5; 
end

% First perform connected component clean up
[ci, sizes, label_id] = CBIG_ComputeConnectedComponentsFromSurface(avg_mesh, labels);
old_labels = labels;
for i = 1:length(label_id)
    conn_sizes = sizes{i};
    for j = 1:length(conn_sizes)
        if(conn_sizes(j) < abs_threshold)
            tmp_index = find(old_labels == label_id(i));
            index = tmp_index(ci{i} == j);
            
            
            neighbors = avg_mesh.vertexNbors(:, index);
            neighbors = reshape(neighbors, numel(neighbors), 1);
            neighbors = neighbors(neighbors ~= 0);
            
            neighbor_labels = old_labels(neighbors); 
            neighbor_labels = neighbor_labels(neighbor_labels ~= label_id(i));
            labels(index) = mode(neighbor_labels(neighbor_labels~=0));
        end
    end
end



function val2 = mode_filter(vertexNbors, val)

temp = vertexNbors;
vertexNbors(vertexNbors == 0) = 1;
val2 = val(vertexNbors);
val2(temp == 0) = min(val) - 1; %ensures no neighbor does not dominate
val2 = mode(val2, 1);

