function [vertices, label, colortable] = CBIG_MapSingleHemiSurfaceDataToAnnotation(data, hemi, ref_annot_dir, discretization_res, colorscale, min_thresh, max_thresh)  
  % [vertices, label, colortable] = CBIG_MapSingleHemiSurfaceDataToAnnotation(data, hemi, ref_annot_dir, discretization_res, colorscale, min_thresh, max_thresh)  
  %
  % Map surface data of one single hemisphere to a format suitable for surface annotation
  %
  % Input:
  %   - data              : original surface data
  %   - hemi              : either 'lh' or 'rh' corresponding to the left or right hemisphere
  %   - ref_annot_dir     : directory containing the reference annotation
  %   - discretization_res: number of distinct discrete values in the output
  %   - colorscale        : the color scheme used for displaying the output data
  %   - min_thresh        : minimum threshold of the output values
  %   - max_thresh        : maximum threshold of the output values
  % Output:
  %   - vertices          : vertices in the reference annotation
  %   - label             : labels assigned to the vertices
  %   - colortable        : colors assigned to the labels 
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  ref_annot_path = fullfile(ref_annot_dir, [hemi '.aparc.annot']);
  [vertices, ref_label, ~] = read_annotation(ref_annot_path);
  
  % set all values above the threshold's upper limit to the threshold's upper limit
  if max_thresh < max(data)
    data(data >= max_thresh) = max_thresh;
  end
  
  % discretize non-zero values
  nonzero_values = data(data >= min_thresh);
  value_step = (max_thresh - min_thresh) / (discretization_res-1);
  binranges = min_thresh:value_step:max_thresh;
  [~, bin_ind] = histc(nonzero_values, binranges);
   
  colortable = CBIG_GenerateAnnotationColortable(discretization_res, colorscale);
  
  underlay_rgb_sum = colortable.table(end, 5);
  
  label = underlay_rgb_sum * ones(length(ref_label), 1); % fill underlay color
  label(data >= min_thresh) = colortable.table(bin_ind, 5);
