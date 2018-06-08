function CBIG_GenerateYeo2011NetworkBoundries
% This function will generate annot files of Yeo2011 7 and 17 network boundaries.
% It can be used as underlay of your surface data.
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
FREESURFER_HOME = getenv('FREESURFER_HOME');

boundary_color = [38 38 38];
boundary_label = boundary_color(1) + boundary_color(2)*2^8 + boundary_color(3)*2^16;

for network_names = {'Yeo2011_7Networks_N1000', 'Yeo2011_17Networks_N1000'}
    network_name = network_names{1};
    for hemis = {'lh' 'rh'}
        hemi = hemis{1};

        % Read fsaverage mesh and Yeo2011 network labels
        mesh = CBIG_ReadNCAvgMesh(hemi, 'fsaverage', 'inflated', [network_name '.annot']);
        labels = mesh.MARS_label;
        
        % Find the boundary if the neighbor verteices have different labels
        BoundaryVec = zeros(length(labels), 1);
        maxNeighbors = size(mesh.vertexNbors, 1);
        for i = 1:length(labels)
            label_vertex = int32(labels(i));
            for k = 1:maxNeighbors
                v_neighbor = mesh.vertexNbors(k, i);
                if(v_neighbor ~= 0 && int32(labels(v_neighbor)) ~= label_vertex)
                    BoundaryVec(i) = 1;
                end
            end
        end
        
        % Write the boundary vector into label and output it into annot file
        [tmp_vertices, tmp_label, tmp_ct] = read_annotation([FREESURFER_HOME '/subjects/fsaverage/label/' hemi '.' network_name '.annot']);
        filename = [hemi '.' network_name '_Boundary.annot'];
        vertices = tmp_vertices;
        label = BoundaryVec * boundary_label;
        ct = tmp_ct;
        ct.numEntries = 1;
        ct.struct_names = {'boundary'};
        ct.table = [boundary_color 0 boundary_label];
        write_annotation(filename, vertices, label, ct);
    end
end
