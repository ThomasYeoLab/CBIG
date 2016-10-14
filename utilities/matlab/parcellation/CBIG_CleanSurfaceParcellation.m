function labels = CBIG_CleanSurfaceParcellation(avg_mesh, labels, num_smooth, connect_frac)

% labels = CBIG_CleanSurfaceParcellation(avg_mesh, labels, num_smooth, connect_frac)
%
% First perform connected component clean up.
%     - The components whose size is smaller that the threshold given by
%     connect_frac will be assigned to the mode label of their neighboring
%     vertices.
%     - If connect_frac is smaller than or equal to 0, this step is skiped.
%
% Second smooth parcellation by mode filter.
%     - The number of iterated smoothing step is given by num_smooth.
%     - If num_smooth is 0, smoothing is skiped.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(nargin < 3)
   num_smooth = 5; 
end

if(nargin < 4)
   connect_frac = 0.1; 
end
abs_threshold = connect_frac*length(labels)/max(labels);

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
            labels(index) = mode(neighbor_labels);
        end
    end
end

for i = 1:num_smooth
  labels = mode_filter(avg_mesh.vertexNbors, labels);
end



function val2 = mode_filter(vertexNbors, val)

temp = vertexNbors;
vertexNbors(vertexNbors == 0) = 1;
val2 = val(vertexNbors);
val2(temp == 0) = min(val) - 1; %ensures no neighbor does not dominate
val2 = mode(val2, 1);

