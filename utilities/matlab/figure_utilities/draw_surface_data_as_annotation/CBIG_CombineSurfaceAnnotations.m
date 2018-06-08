function CBIG_CombineSurfaceAnnotations(underlay_annot_file, overlay_annot_file, combined_annot_file, excluded_label)

% CBIG_CombineSurfaceAnnotations(underlay_annot_file, overlay_annot_file, combined_annot_file, excluded_label)
%
% Combine two surface annotations
%
% Input:
%   Compulsory:
%     - underlay_annot_file: file containing the underlaying annotation
%     - overlay_annot_file : file containing the overlaying annotation
%     - combined_annot_file: file containing the combined annotations
%   Optional:
%     - excluded_label: vertices of the underlaying annotation with label
%     equals to excluded_label would not be overwritten.
%
% Example:
%  CBIG_CombineSurfaceAnnotations('/data/Work/lh_component.annot', '/data/Work/lh_parcel_outlines.annot', ...
%              '/data/Work/lh_component_with_parcel_outlines.annot', 0);
%  Overlay the surface annotation saved in
%  '/data/Work/lh_parcel_outlines.annot' on top of the surface annotation
%  saved in '/data/Work/lh_component.annot' and output the new surface
%  annotation in '/data/Work/lh_component_with_parcel_outlines.annot'. Vertices
%  with label 0 in '/data/Work/lh_component.annot' are  not overwritten.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    [vertices, first_labels, first_colortable] = read_annotation(underlay_annot_file);
    [tmp, second_labels, second_colortable] = read_annotation(overlay_annot_file);
    
    combined_vertices = vertices;
    combined_labels = first_labels;
    combined_colortable = first_colortable;
    
    numEntries = combined_colortable.numEntries;
    for idx = second_colortable.numEntries
        color = second_colortable.table(idx, :);
        colorcode = color(5);
        
        vertices_of_curr_colorcode = (second_labels == colorcode);
        % Exclude vertices with label = excluded_label
        if nargin == 4
          vertices_of_curr_colorcode(first_labels == excluded_label) =0;
        end
        
        numEntries = numEntries + 1;
        combined_colortable.numEntries = numEntries;
        combined_colortable.struct_names{numEntries} = second_colortable.struct_names{idx};
        
        combined_labels(vertices_of_curr_colorcode) = colorcode;
        combined_colortable.table(numEntries, :) = color;
        if (sum(combined_colortable.table(1:numEntries-1, 5) == colorcode) > 0)
            max_colorcode = max(combined_colortable.table(:, 5));
            new_colorcode = max_colorcode + 1;
            
            combined_colortable.table(numEntries, 5) = new_colorcode;
            combined_labels(vertices_of_curr_colorcode) = new_colorcode;
        end
    end
    
    write_annotation(combined_annot_file, combined_vertices, combined_labels, combined_colortable);
