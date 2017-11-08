function surface_normals = CBIG_ComputeNormalVectorsCaret(mesh, avg_neighbor_bool)

% surface_normals = CBIG_ComputeNormalVectorsCaret(mesh, avg_neighbor_bool)
% This function is used to compute normal vectors for each vertex in mesh.
% surface_normals = 3 x # vertices 
% assumes input mesh is a "regular" subdivided icosahedron
% # vertices = V, # faces = F
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(nargin < 2)
    avg_neighbor_bool = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute unit normals for each triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vertexAcoords = mesh.vertices(:, mesh.faces(1, :));
vertexBcoords = mesh.vertices(:, mesh.faces(2, :));
vertexCcoords = mesh.vertices(:, mesh.faces(3, :));
vecAB = vertexBcoords - vertexAcoords;
vecAC = vertexCcoords - vertexAcoords;

% "minus" forces normals to point outwards.
% to test, read spherical mesh
% tmp = dot(mesh.vertices(:,  mesh.faces(1, :)), triangle_normals, 1);
% min(tmp) should be greater than 0.
triangle_normals = -cross(vecAB, vecAC, 1); % 3 x F
triangle_normals = bsxfun(@times, triangle_normals, 1./sqrt(sum(triangle_normals.^2, 1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute unit normals for each vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valid_neighbor_triangles = (mesh.vertexFaces ~= 0);
num_valid_neighbor_triangles = sum(valid_neighbor_triangles, 1);
neighbor_triangles = mesh.vertexFaces; neighbor_triangles(neighbor_triangles == 0) = 1;

normal_x = triangle_normals(1, neighbor_triangles); 
normal_x(~valid_neighbor_triangles) = 0;
normal_x = reshape(normal_x, size(neighbor_triangles));
surface_normals(1, :) = sum(normal_x, 1)./num_valid_neighbor_triangles;

normal_y = triangle_normals(2, neighbor_triangles); 
normal_y(~valid_neighbor_triangles) = 0;
normal_y = reshape(normal_y, size(neighbor_triangles));
surface_normals(2, :) = sum(normal_y, 1)./num_valid_neighbor_triangles;

normal_z = triangle_normals(3, neighbor_triangles); 
normal_z(~valid_neighbor_triangles) = 0;
normal_z = reshape(normal_z, size(neighbor_triangles));
surface_normals(3, :) = sum(normal_z, 1)./num_valid_neighbor_triangles;

surface_normals = bsxfun(@times, surface_normals, 1./sqrt(sum(surface_normals.^2, 1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth normal if avg_neighbor_bool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(avg_neighbor_bool)
    surface_normals = MARS_AverageData(mesh, single(surface_normals), 0, 1);
    surface_normals = bsxfun(@times, surface_normals, 1./sqrt(sum(surface_normals.^2, 1)));
end

