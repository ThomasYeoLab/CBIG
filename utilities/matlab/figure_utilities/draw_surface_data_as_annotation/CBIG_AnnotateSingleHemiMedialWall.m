function CBIG_AnnotateSingleHemiMedialWall(vertices, labels, colortable, ref_labels, ref_medialwall_label, output_annot_file)

% CBIG_AnnotateSingleHemiMedialWall(vertices, labels, colortable, ref_labels, ref_medialwall_label, output_annot_file)
%
% Mark the medial wall with dark gray color on a surface annotation on
% a single hemisphere. The medial wall is delineated in the reference
% annotation file. The resulting annotation is written to an output file.
%
% Input:
%   - vertices            : vertices of the annotation
%   - labels              : original labels of the annotation
%   - colortable          : original colortable of the annotation
%   - ref_labels          : labels of the reference annotation
%   - ref_medialwall_label: label assigned to the medial wall in the reference
%                           annotation
%   - output_annot_file   : output surface annotation 
%
% Example:
%   abs_path_to_ref_lh_annot_file = fullfile(getenv('FREESURFER_HOME'), ...
%                  'subjects', 'fsaverage5', 'label', 'lh.aparc.annot');
%   [lh_vertices, ref_lh_labels, ~] = read_annotation(abs_path_to_ref_lh_annot_file);
%   ref_medial_wall_albel = 0;
%   CBIG_AnnotateSingleHemiMedialWall(lh_vertices, lh_labels, lh_colortable, ...
%          ref_lh_labels, 0, '/data/Work/lh_with_medial_wall.annot');
%
%   Draw a meidal wall ontop of the original annotation saved in lh_labels.
%   The colortable of the original annotation is saved in lh_colortable.
%   The reference vertices and reference labels are loaded in lh_vertices
%   and ref_lh_labels respetively (which comes from the default fsaverage5
%   aparc.annot parcellation in this example. The medialwall is labelled as
%   0 in the reference annotation. The output annotation is saved at
%   '/data/Work/lh_with_medial_wall.annot'
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    % pick out the medial wall and define its new color
    medialwall_vertices = (ref_labels == ref_medialwall_label);
    medialwall_color = zeros(1, 5);
    medialwall_color(1:3) = [50 50 50];
    medialwall_color(5)   = 50 + 50 * 2^8 + 50 *2^16;
    
    labels(medialwall_vertices) = medialwall_color(5);
    
    % append the medial wall's color to the colortable of the annotation
    new_colortable_table = colortable;
    if ~any(new_colortable_table.table(:, 5) == medialwall_color(5))
        newNumEntries = colortable.numEntries + 1;
        new_colortable_table.numEntries = newNumEntries;
        new_colortable_table.struct_names{newNumEntries} = 'medialwall';
        new_colortable_table.table(newNumEntries, :) = medialwall_color;
    end
    
    write_annotation(output_annot_file, vertices, labels, new_colortable_table);
