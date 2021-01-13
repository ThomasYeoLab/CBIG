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
function MARS_sbjMesh = MARS2_readSbjMesh(params)

% params:
% {SUBJECTS_DIR, SUBJECT, surf_filename, metric_surf_filename, 
% annot_filename, data_filename_cell, unfoldBool, structures_of_interest,
% normalizeBool, hemi, radius}

%This function read in a surface specified by surf_filename found in the
%SUBJECTS_DIR/SUBJECT/surf directory, as well as the annotation by
%annot_filename found in the SUBJECTS_DIR/SUBJECT/label directory and
%finally the data specified by the data_filename_cell assumed to be in the
%SUBJECTS_DIR/SUBJECT/label directory.
%function also reads in surface for which metric distortion is calculated
%specified by metric_surf_filename
%
%Assume normalization at 100.

if(~isfield(params, 'radius'))
    radius = 100;
else
    radius = params.radius;
end


total_surface_area = 4*pi*(radius^2);

% Read surface specified by user
if(isfield(params, 'hemi'))
    [vertices, faces] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.hemi '.' params.surf_filename]);
else
    [vertices, faces] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.surf_filename]);
end

%Reverse direction of faces!!! so as to be in line with icosahedronMesh
if (isfield(params, 'flipFacesBool'))
    if(params.flipFacesBool)
        faces = [faces(:,3) faces(:,2) faces(:,1)];
    end
end

vertices = single(vertices);
faces = int32(faces+1);

LengthVec = sqrt(sum(vertices.^2, 2));
vertices = vertices./repmat(LengthVec, 1, 3) * radius; %normalizes surface area to that specified by radius

% Read the surface used to calculate metric distortion

if (isfield(params, 'metric_surf_filename'))
    if(isfield(params, 'hemi'))
        [metricVerts, ignored] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.hemi '.' params.metric_surf_filename]);
    else
        [metricVerts, ignored] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.metric_surf_filename]);
    end
    metricVerts = single(metricVerts);
    %normalize surface area to be the same.
    %metricSurfaceArea = MARS_calculateSurfaceArea(metricVerts, faces);
    
    metricSurfaceArea = sum(MARS_computeMeshFaceAreas(int32(size(faces, 1)), int32(faces'), single(metricVerts')), 'double');
    metricSurfaceArea = single(metricSurfaceArea);
    surface_scaling_factor = sqrt(total_surface_area/metricSurfaceArea);
    metricVerts = metricVerts*surface_scaling_factor;
else
    surface_scaling_factor = [];
    metricVerts = [];
end


%Read data
if (isfield(params, 'data_filename_cell'))
    num_data = length(params.data_filename_cell);
    data = zeros(length(vertices), num_data, 'single');

    normalization = zeros(num_data, 1);
    for i = 1:num_data
        if(isfield(params, 'hemi'))
            data(:,i) = single(read_curv([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.hemi '.' params.data_filename_cell{i}]));
        else
            data(:,i) = single(read_curv([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.data_filename_cell{i}]));
        end
        if(isfield(params, 'normalizeBool'))
            if(params.normalizeBool == 1)
                data(:,i) = data(:,i)-median(data(:,i));
                data(:,i) = data(:,i)/std(data(:,i));
                
                index = find(data(:, i) < -3);
                data(index, i) = -3 - (1 - exp(3 - abs(data(index, i))));
                index = find(data(:, i) > 3);
                data(index, i) = 3 + (1 - exp(3 - abs(data(index, i))));
                
                data(:,i) = data(:,i)/std(data(:,i));
                index = find(data(:, i) < -3);
                data(index, i) = -3 - (1 - exp(3 - abs(data(index, i))));
                index = find(data(:, i) > 3);
                data(index, i) = 3 + (1 - exp(3 - abs(data(index, i))));
            elseif(params.normalizeBool == 2)   
                
                %max normalization
                BIG_NO = -1e6;
                if(data(1, i) ~= BIG_NO)
                    normalization(i) = max(data(:,i));
                    data(:, i) = data(:, i)/normalization(i);
                end
                
            else
                clear normalization;%do nothing
            end
        end
    end
else
    data = [];
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
if (isfield(params, 'annot_filename'))
    if(isfield(params, 'hemi'))
        [ignored, MARS_label, MARS_ct] = read_annotation([params.SUBJECTS_DIR '/' params.SUBJECT '/label/' params.hemi '.' params.annot_filename]);
    else
        [ignored, MARS_label, MARS_ct] = read_annotation([params.SUBJECTS_DIR '/' params.SUBJECT '/label/' params.annot_filename]);
    end

    if (isfield(params, 'structures_of_interest'))
        MARS_ct = MARS_reorganizeCT(MARS_ct, params.structures_of_interest);
        
    end
    [MARS_label, MARS_ct] = MARS_reorganizeLabels(MARS_label, MARS_ct, vertexNbors');
else
    MARS_label = [];
    MARS_ct = [];    
end

% Read Brain labels
if(isfield(params, 'label_filename_cell'))
   if(isfield(params, 'annot_filename'))
     error('Should not have both annotation files and label files'); 
   end
    
   MARS_ct.orig_tab = 'FunctionalCoordinate';
   MARS_ct.numEntries = 1;
   MARS_ct.struct_names = {'Unknown'};
   MARS_ct.table = [0 0 0 0 0 1];
   
   MARS_label = ones(size(vertices, 1), 1);
   for i = 1:length(params.label_filename_cell)
       l = read_label([], fullfile(params.SUBJECTS_DIR, params.SUBJECT, 'label', [params.hemi '.' params.label_filename_cell{i} '.label']));
       MARS_label(l(:, 1)+1) = i+1;
       
       MARS_ct.numEntries = MARS_ct.numEntries + 1;
       MARS_ct.struct_names{i+1} = params.label_filename_cell{i};
       MARS_ct.table = [MARS_ct.table; ...
                        i i i 0 i+i*2^8+i*2^16 i+1];
   end
   MARS_label = int32(MARS_label);
   MARS_ct.table = int32(MARS_ct.table);
   MARS_ct.numEntries = int32(MARS_ct.numEntries);
end


%Find vertexDistSq2Nbors
vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(vertexNbors', 1)), int32(size(vertices', 2)), int32(vertexNbors'), single(vertices'));

% Note the transpose!!!!!
MARS_sbjMesh = struct('vertices', vertices', 'metricVerts', metricVerts', 'faces', faces', 'vertexNbors', vertexNbors', 'vertexFaces', vertexFaces', 'vertexDistSq2Nbors', vertexDistSq2Nbors, ...
    'faceAreas', faceAreas', 'data', data', 'MARS_label', MARS_label', 'MARS_ct', MARS_ct, 'surface_scaling_factor', surface_scaling_factor);
if (isfield(params, 'unfoldBool'))
    if(params.unfoldBool)        
        if (isfield(params, 'flipFacesBool'))
            MARS_sbjMesh = MARS_simpleUnfoldMesh(MARS_sbjMesh, ~params.flipFacesBool);
        else
            MARS_sbjMesh = MARS_simpleUnfoldMesh(MARS_sbjMesh, 1); %this is for backward compatibility.
        end
    end
end

% Read second surface if necessary
if(isfield(params, 'surf_filename2'))
    if(isfield(params, 'hemi'))
        [vertices2] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.hemi '.' params.surf_filename2]);
    else
        [vertices2] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.surf_filename2]);
    end
    
    vertices2 = single(vertices2);
    LengthVec = sqrt(sum(vertices2.^2, 2));
    vertices2 = vertices2./repmat(LengthVec, 1, 3) * radius; %normalizes surface area to that specified by radius
    vertices2 = vertices2';

    if(size(vertices2, 2) ~= size(MARS_sbjMesh.vertices, 2))
       error('Num of vertices in surf_filename2 not the same as surf_filename'); 
    end
    
    temp_sbjMesh = MARS_sbjMesh;
    temp_sbjMesh.vertices = vertices2;
    if(isfield(params, 'unfoldBool'))
        if(params.unfoldBool)
            if (isfield(params, 'flipFacesBool'))
                temp_sbjMesh = MARS_simpleUnfoldMesh(temp_sbjMesh, ~params.flipFacesBool);
            else
                temp_sbjMesh = MARS_simpleUnfoldMesh(temp_sbjMesh, 1); %this is for backward compatibility.
            end
        end
    end
    MARS_sbjMesh.vertices2 = temp_sbjMesh.vertices;
    
    MARS_sbjMesh.vertexDistSq2Nbors2 = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)), int32(size(MARS_sbjMesh.vertices2, 2)), ...
                                                                      int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.vertices2));   
end

% Read second set of data_filename_cell2
if (isfield(params, 'data_filename_cell2'))
    num_data = length(params.data_filename_cell2);
    data = zeros(length(vertices), num_data, 'single');

    normalization = zeros(num_data, 1);
    for i = 1:num_data
        if(isfield(params, 'hemi'))
            data(:,i) = single(read_curv([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.hemi '.' params.data_filename_cell2{i}]));
        else
            data(:,i) = single(read_curv([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.data_filename_cell2{i}]));
        end
        if(isfield(params, 'normalizeBool2'))
            if(params.normalizeBool2 == 1)
                data(:,i) = data(:,i)-median(data(:,i));
                data(:,i) = data(:,i)/std(data(:,i));
                
                index = find(data(:, i) < -3);
                data(index, i) = -3 - (1 - exp(3 - abs(data(index, i))));
                index = find(data(:, i) > 3);
                data(index, i) = 3 + (1 - exp(3 - abs(data(index, i))));
                
                data(:,i) = data(:,i)/std(data(:,i));
                index = find(data(:, i) < -3);
                data(index, i) = -3 - (1 - exp(3 - abs(data(index, i))));
                index = find(data(:, i) > 3);
                data(index, i) = 3 + (1 - exp(3 - abs(data(index, i))));
            elseif(params.normalizeBool2 == 2)   
                
                %max normalization
                BIG_NO = -1e6;
                if(data(1, i) ~= BIG_NO)
                    normalization(i) = max(data(:,i));
                    data(:, i) = data(:, i)/normalization(i);
                end
            else
                clear normalization;%do nothing
            end
        end
    end
    MARS_sbjMesh.data2 = transpose(data);
    
    if(size(MARS_sbjMesh.data2, 1) == 1) %if only read in one DT, cheat to produce known DT
        MARS_sbjMesh.data2 = [-MARS_sbjMesh.data2; MARS_sbjMesh.data2];
    end
end

% Read third surface if necessary
if(isfield(params, 'surf_filename3'))
    if(isfield(params, 'hemi'))
        [vertices3] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.hemi '.' params.surf_filename3]);
    else
        [vertices3] = read_surf([params.SUBJECTS_DIR '/' params.SUBJECT '/surf/' params.surf_filename3]);
    end
    
    vertices3 = single(vertices3);
    LengthVec = sqrt(sum(vertices3.^2, 2));
    vertices3 = vertices3./repmat(LengthVec, 1, 3) * radius; %normalizes surface area to that specified by radius
    vertices3 = vertices3';

    if(size(vertices3, 2) ~= size(MARS_sbjMesh.vertices, 2))
       error('Num of vertices in surf_filename3 not the same as surf_filename'); 
    end
    
    temp_sbjMesh = MARS_sbjMesh;
    temp_sbjMesh.vertices = vertices3;
    if(isfield(params, 'unfoldBool'))
        if(params.unfoldBool)
            if (isfield(params, 'flipFacesBool'))
                temp_sbjMesh = MARS_simpleUnfoldMesh(temp_sbjMesh, ~params.flipFacesBool);
            else
                temp_sbjMesh = MARS_simpleUnfoldMesh(temp_sbjMesh, 1); %this is for backward compatibility.
            end
        end
    end
    MARS_sbjMesh.vertices3 = temp_sbjMesh.vertices;
    
    MARS_sbjMesh.vertexDistSq2Nbors3 = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)), int32(size(MARS_sbjMesh.vertices3, 2)), ...
                                                                      int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.vertices3));   
end

% Create stucture existence matrix
if(~isempty(MARS_sbjMesh.MARS_label))
    MARS_sbjMesh.struct_exist = zeros(MARS_sbjMesh.MARS_ct.numEntries, 1);
    for i = 1:length(MARS_sbjMesh.struct_exist)
        if(~isempty(find(MARS_sbjMesh.MARS_label==i, 1)))
            MARS_sbjMesh.struct_exist(i) = 1;
        end
    end
end

