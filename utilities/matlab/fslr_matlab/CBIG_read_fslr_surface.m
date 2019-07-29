function output_mesh = CBIG_read_fslr_surface(hemi, mesh_name, surf_type, label)

% output_mesh = CBIG_read_fslr_surface(hemi, mesh_name, surf_type, label)
%
% This is an fs_LR version for CBIG_ReadNCAvgMesh/MARS2_readSbjMesh 
% Input:
%      - hemi:
%        'lh' or 'rh' for left or right hemisphere respectively
%      - mesh_name:
%        mesh structure. mesh_name = 'fs_LR_32k' or 'fs_LR_164k'.
%      - surf_type: 
%        surface template. Options are 'inflated', 'very_inflated', 
%        'midthickness_orig', 'white_orig', 'pial_orig', 'sphere'.
%      - label: 
%        If label is not defined, then there is no 'MARS_label' and
%        'MARS_ct' structure fields in the output structure. label can ONLY
%        be set to 'aparc.annot' or 'medialwall.annot'.
%        'aparc.annot' is the Desikan parcellation in fslr.
%        'medialwall.annot' is used to look up medial wall vertices 
%        (output_mesh.MARS_label will be 1 for medial wall vertices, and 2 for cortical vertices.)
% Output:
%      - output_mesh:
%        output mesh structure.
%
% Example:
% lh_mesh_fslr_164k = CBIG_read_fslr_surface('lh','fs_LR_164k','inflated','aparc.annot');
% lh_mesh_fslr_32k  = CBIG_read_fslr_surface('lh','fs_LR_32k','very_inflated');
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


    if(~exist('mesh_name', 'var'))
       mesh_name = 'fs_LR_32k'; 
    end
    if(~exist('surf_type', 'var'))
       surf_type = 'inflated'; 
    end
    if(nargin == 4)
        if(strfind(label, 'annot'))
            annot_filename = label;
        end
    end
    
    
    % Read surface specified by user
    folder = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface',mesh_name);
    radius = 100;
    total_surface_area = 4*pi*(radius^2);
    [vertices,faces] = read_surf(fullfile(folder, 'surf', [hemi, '.', surf_type]));
    output_mesh.faces = int32(faces')+1;
    output_mesh.vertices = single(vertices');
    
    
    %%% compute vertexNbors
    output_mesh.vertexNbors = MARS_convertFaces2VertNbors(output_mesh.faces', ...
        int32(size(output_mesh.vertices,2)));
    num_per_vertex = length(output_mesh.vertexNbors)/size(output_mesh.vertices,2);
    output_mesh.vertexNbors = reshape(output_mesh.vertexNbors, ...
        size(output_mesh.vertices,2), num_per_vertex);
    output_mesh.vertexNbors=output_mesh.vertexNbors';


    %%% compute metricVerts
    metricSurfaceArea = sum(MARS_computeMeshFaceAreas(int32(size(output_mesh.faces, 2)), ...
        int32(output_mesh.faces), single(output_mesh.vertices)));
    output_mesh.surface_scaling_factor = sqrt(total_surface_area/metricSurfaceArea);
    output_mesh.metricVerts = output_mesh.vertices*output_mesh.surface_scaling_factor;
    
    
    %Find vertexFaces
    output_mesh.vertexFaces =  MARS_convertFaces2FacesOfVert(output_mesh.faces', ...
        int32(size(output_mesh.vertices, 2)));
    output_mesh.vertexFaces = reshape(output_mesh.vertexFaces, ...
        size(output_mesh.vertices,2), num_per_vertex);
    output_mesh.vertexFaces = output_mesh.vertexFaces';
    
    
    % Compute Face Areas.
    output_mesh.faceAreas = MARS_computeMeshFaceAreas(int32(size(output_mesh.faces, 2)), ...
        int32(output_mesh.faces), single(output_mesh.vertices));
    output_mesh.faceAreas = output_mesh.faceAreas';
    
    
    %Read brain labels
    if (exist('annot_filename', 'var'))
        [~, MARS_label, MARS_ct] = read_annotation(fullfile(folder, ...
            'label', [hemi, '.', annot_filename]));
        [output_mesh.MARS_label, output_mesh.MARS_ct] = ...
            MARS_reorganizeLabels(MARS_label, MARS_ct, output_mesh.vertexNbors');  
    end

end
