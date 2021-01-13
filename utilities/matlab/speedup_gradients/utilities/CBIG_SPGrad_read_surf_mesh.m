function sbjMesh = CBIG_SPGrad_read_surf_mesh(surf_path)

% This function reads a surface mesh from a surface file
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


radius = 100;
total_surface_area = 4*pi*(radius^2);
[vertices, faces] = read_surf(surf_path);

vertices = single(vertices);
faces = int32(faces+1);

metricVerts = single(vertices);
metricSurfaceArea = sum(MARS_computeMeshFaceAreas(int32(size(faces, 1)), int32(faces'), single(metricVerts')));

surface_scaling_factor = sqrt(total_surface_area/metricSurfaceArea);
metricVerts = metricVerts*surface_scaling_factor;
vertexFaces =  MARS_convertFaces2FacesOfVert(faces, int32(size(vertices, 1)));
num_per_vertex = length(vertexFaces)/size(vertices,1);

vertexFaces = reshape(vertexFaces, size(vertices,1), num_per_vertex);
faceAreas = MARS_computeMeshFaceAreas(int32(size(faces, 1)), int32(faces'), single(vertices'));  
vertexNbors = MARS_convertFaces2VertNbors(faces, int32(size(vertices,1)));
num_per_vertex = length(vertexNbors)/size(vertices,1);

vertexNbors = reshape(vertexNbors, size(vertices,1), num_per_vertex);
vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(vertexNbors', 1)), ...
int32(size(vertices', 2)), int32(vertexNbors'), single(vertices'));

data = [];MARS_label = []; MARS_ct = []; 

sbjMesh = struct('vertices', vertices', 'metricVerts', metricVerts', 'faces', faces', ...
    'vertexNbors', vertexNbors', 'vertexFaces', vertexFaces', 'vertexDistSq2Nbors', ...
    vertexDistSq2Nbors, 'faceAreas', faceAreas', 'data', data', ...
    'MARS_label', MARS_label', 'MARS_ct', MARS_ct, 'surface_scaling_factor', surface_scaling_factor);
end