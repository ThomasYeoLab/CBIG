function CBIG_AnnotateMedialWall(annot_file, hemi, output_annot_file)

  % CBIG_AnnotateMedialWall(annot_file, hemi, output_annot_file)
  %
  % Mark the medial wall with dark gray color on a surface annotation
  %
  % Input:
  %   - annot_file       : original surface annotation file
  %   - hemi             : either 'lh' or 'rh' corresponding to the left and right hemisphere
  %   - output_annot_file: output surface annottion 
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  [vertices, labels, colortable] = read_annotation(annot_file);

  if length(vertices) == 40962 % fsaverage6
    ref_annot_dir = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage6', 'label');
  elseif length(vertices) == 10242 % fsaverage5
    ref_annot_dir = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage5', 'label');
  elseif length(vertices) == 163842 % fsaverage
    ref_annot_dir = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage', 'label');
  end
  
  ref_annot_path = fullfile(ref_annot_dir, [hemi '.aparc.annot']);
  [~, ref_label, ~] = read_annotation(ref_annot_path);
  medialwall_vertices = (ref_label == 0);
  medialwall_color = zeros(1, 5);
  medialwall_color(1:3) = [50 50 50];
  medialwall_color(5)   = 50 + 50 * 2^8 + 50 *2^16;
  
  new_table = colortable;
  if ~any(new_table.table(:, 5) == medialwall_color(5))
    newNumEntries = colortable.numEntries + 1;
    new_table.numEntries = newNumEntries;
    new_table.struct_names{newNumEntries} = 'medialwall';
    new_table.table(newNumEntries, :) = medialwall_color;
  end
  
  labels(medialwall_vertices) = medialwall_color(5);
  write_annotation(output_annot_file, vertices, labels, new_table);
