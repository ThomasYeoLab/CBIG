function neighborhood_cell=CBIG_ComputeNeighborhood(hemi, mesh_file, radius, surf_type)

% neighborhood_cell=CBIG_ComputeNeighborhood(hemi, mesh_file, radius, surf_type)
% This function is used to find neighborhood vertices for each vertex
% within radius.
% Input:
%      -hemi: 'lh' or 'rh'
%      -mesh_file: Freesurfer mesh structure. e.g. 'fsaverage5'
%      -radius: geometric distance defined by shortest path
%      -surf_type: 'inflated', 'sphere', or 'white'
% Example:
% CBIG_ComputeNeighborhood('lh', 'fsaverage5', '12', 'sphere'); CBIG_ComputeNeighborhood('rh', 'fsaverage5', '12', 'sphere');
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(radius))
   radius = str2num(radius); 
end

avg_mesh = CBIG_ReadNCAvgMesh(hemi, mesh_file, surf_type);
n = size(avg_mesh.vertices, 2);
g = zeros(n, n);
for i = 1:n
    for j = 1:size(avg_mesh.vertexNbors, 1) % foreach neighbor
        if(avg_mesh.vertexNbors(j, i) > 0)
            g(i, avg_mesh.vertexNbors(j, i)) = sqrt(avg_mesh.vertexDistSq2Nbors(j, i));
        end
    end
end
g = sparse(g);
D_geo  = all_shortest_paths(g);

neighborhood_cell = cell(n, 1);
for i = 1:n
   neighborhood_cell{i} = find(D_geo(i, :) <= radius); 
end
save([hemi '.' mesh_file '.' num2str(radius) '.' surf_type '.mat'], 'neighborhood_cell');



