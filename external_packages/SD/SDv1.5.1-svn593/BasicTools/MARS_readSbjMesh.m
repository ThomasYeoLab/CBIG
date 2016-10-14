% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================
function MARS_sbjMesh = MARS_readSbjMesh(SUBJECTS_DIR, SUBJECT, surf_filename, metric_surf_filename, annot_filename, data_filename_cell, unfoldBool, structures_of_interest, normalizeBool, radius)

% MARS_sbjMesh = MARS_readSbjMesh(SUBJECTS_DIR, SUBJECT, surf_filename, metric_surf_filename, annot_filename, data_filename_cell, unfoldBool, structures_of_interest, normalizeBool, radius)
%
%This function read in a surface specified by surf_filename found in the
%SUBJECTS_DIR/SUBJECT/surf directory, as well as the annotation by
%annot_filename found in the SUBJECTS_DIR/SUBJECT/label directory and
%finally the data specified by the data_filename_cell assumed to be in the
%SUBJECTS_DIR/SUBJECT/label directory.
%function also reads in surface for which metric distortion is calculated
%specified by metric_surf_filename
%
%Assume normalization at 100.

if(nargin < 10)
    radius = 100;
end
total_surface_area = 4*pi*radius^2;

% Read surface specified by user
[vertices, faces] = read_surf([SUBJECTS_DIR '/' SUBJECT '/surf/' surf_filename]);

%Reverse direction of faces!!! so as to be in line with icosahedronMesh
%faces = [faces(:,3) faces(:,2) faces(:,1)];

vertices = single(vertices);
faces = int32(faces+1);

LengthVec = sqrt(sum(vertices.^2, 2));
vertices = vertices./repmat(LengthVec, 1, 3) * radius; %normalizes surface area to that specified by radius

% Read the surface used to calculate metric distortion
[metricVerts, ignored] = read_surf([SUBJECTS_DIR '/' SUBJECT '/surf/' metric_surf_filename]);
metricVerts = single(metricVerts);

%normalize surface area to be the same.
metricSurfaceArea = MARS_calculateSurfaceArea(metricVerts, faces);
surface_scaling_factor = sqrt(total_surface_area/metricSurfaceArea);
metricVerts = metricVerts*surface_scaling_factor;

%Read data
num_data = length(data_filename_cell);
data = zeros(length(vertices), num_data, 'single');

for i = 1:num_data
    data(:,i) = single(read_curv([SUBJECTS_DIR '/' SUBJECT '/surf/' data_filename_cell{i}]));
    if(normalizeBool)
        data(:,i) = data(:,i)-mean(data(:,i));
        data(:,i) = data(:,i)/std(data(:,i));
    end
end

%Find vertexFaces
vertexFaces =  MARS_convertFaces2FacesOfVert(faces, int32(size(vertices, 1)));
num_per_vertex = length(vertexFaces)/size(vertices,1);
vertexFaces = reshape(vertexFaces, size(vertices,1), num_per_vertex);

% Compute Face Areas.
faceAreas = MARS_computeMeshFaceAreas(int32(size(faces, 1)), int32(faces'), single(vertices'));  

%Find vertexNbors
vertexNbors = MARS_convertFaces2VertNbors(faces, int32(size(vertices,1)));
num_per_vertex = length(vertexNbors)/size(vertices,1);
vertexNbors = reshape(vertexNbors, size(vertices,1), num_per_vertex);

%Read brain labels
[ignored, MARS_label, MARS_ct] = read_annotation([SUBJECTS_DIR '/' SUBJECT '/label/' annot_filename]);
MARS_ct = MARS_reorganizeCT(MARS_ct, structures_of_interest);
[MARS_label, MARS_ct] = MARS_reorganizeLabels(MARS_label, MARS_ct, vertexNbors');

% Note the transpose!!!!!
MARS_sbjMesh = struct('vertices', vertices', 'metricVerts', metricVerts', 'faces', faces', 'vertexNbors', vertexNbors', 'vertexFaces', vertexFaces', ...
    'faceAreas', faceAreas', 'data', data', 'MARS_label', MARS_label', 'MARS_ct', MARS_ct, 'surface_scaling_factor', surface_scaling_factor);

if(unfoldBool)
    MARS_sbjMesh = MARS_simpleUnfoldMesh(MARS_sbjMesh, 1);
end


