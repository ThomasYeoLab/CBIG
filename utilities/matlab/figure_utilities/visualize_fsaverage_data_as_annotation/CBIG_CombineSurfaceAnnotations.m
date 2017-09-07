function CBIG_CombineSurfaceAnnotations(underlay_annot_file, overlay_annot_file, combined_annot_file)

  % CBIG_CombineSurfaceAnnotations(underlay_annot_file, overlay_annot_file, combined_annot_file)
  %
  % Combine two surface annotations
  %
  % Input:
  %   - underlay_annot_file: file containing the underlaying annotation
  %   - overlay_annot_file : file containing the overlaying annotation
  %   - combined_annot_file: file containing the combined annotations
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
