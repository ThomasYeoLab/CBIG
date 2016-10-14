function avg_mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, surf_type, label)

% Read the average mesh (from FREESURFER) of normal control people.
% 
%   avg_mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, surf_type, label)
%   Input:
%       hemi        : 'lh' or 'rh'
%       mesh_name   : 'fsaverage6', 'fsaverage5', 'fsaverage4'
%       surf_type   : 'inflated', 'white', 'sphere'
%       label       : 'cortex', 'Yeo2011_7Networks_N1000.annot'
%   Output: 
%       avg_mesh    : average mesh
% 
%   PS: For more labels, check $FREESURFER_HOME/subjects/fsaverage/label
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(nargin < 3)
   surf_type = 'white'; 
end

parms.SUBJECTS_DIR = fullfile(getenv('FREESURFER_HOME'), 'subjects');
parms.hemi = hemi;
parms.SUBJECT = mesh_name;
parms.read_surface = @MARS2_readSbjMesh;
parms.radius = 100;
parms.unfoldBool = 0;
parms.flipFacesBool = 1;
parms.surf_filename = 'sphere';
parms.metric_surf_filename = surf_type;
parms.data_filename_cell = {'sulc', 'curv'};

if(nargin == 4)
    if(strfind(label, 'annot'))
        parms.annot_filename = label;
    else
        parms.label_filename_cell = {label};
    end
end

avg_mesh = MARS2_readSbjMesh(parms);
avg_mesh.vertices = avg_mesh.metricVerts;

% Recompute vertexDistSq2Nbors and faceAreas
avg_mesh.vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(avg_mesh.vertexNbors, 1)), int32(size(avg_mesh.vertices, 2)), int32(avg_mesh.vertexNbors), single(avg_mesh.vertices));
avg_mesh.faceAreas = MARS_computeMeshFaceAreas(int32(size(avg_mesh.faces, 2)), int32(avg_mesh.faces), single(avg_mesh.vertices));
