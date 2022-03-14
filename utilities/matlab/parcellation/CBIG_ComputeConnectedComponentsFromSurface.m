function [ci, sizes, label_id] = CBIG_ComputeConnectedComponentsFromSurface(avg_mesh, labels)

% Given average mesh and labels, compute connected components.
% 
%   [ci, sizes, label_id] = CBIG_ComputeConnectedComponentsFromSurface(avg_mesh, labels)
%   Input:
%       avg_mesh    : average surface mesh read from CBIG_ReadNCAvgMesh.m
%       OR CBIG_read_fslr_surface.m
%       labels      : parcellation surface labels
%   Output:
%       ci          : an array of cell of length k, each cell specifies the
%       indices of each connected component
%       sizes       : an array of cell of length k; each cell specifies the sizes of each connected component
%       label_id    : a vector of length k; all unique indices of input labels 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% To iterate through all the components of the parcellation:
%
% for ii = 1:length(label_id)
%   i = label_id(ii);
%   label_vertices = find(labels == i);
%   for j = 1:length(sizes{ii})
%       component_vertices = label_vertices(ci{ii} == j);       
%   end
% end
%

num_labels = max(labels);

% form graph of surface
neighbors = avg_mesh.vertexNbors;
self = repmat(1:length(labels), size(neighbors, 1), 1);

neighbors = reshape(neighbors, numel(neighbors), 1);
self = reshape(self, numel(self), 1);

% remove zero neighbors
non_zero_index = neighbors ~= 0;
neighbors = neighbors(non_zero_index);
self = self(non_zero_index);

g = sparse(double(neighbors), self, ones(size(self)), length(labels), length(labels));

count = 0;
for i = 0:num_labels
    index = find(labels == i);
    if(~isempty(index))
        count = count + 1;
        subgraph = g(index, index);
        [ci{count}, sizes{count}] = components(subgraph);
        label_id(count) = i;
    end
end

