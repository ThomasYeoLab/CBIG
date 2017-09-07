function [lh_annot_file, rh_annot_file] = CBIG_ProjectFsaverageDataToDiscretizedAnnotation(lh_data, rh_data, discretization_res, min_thresh, max_thresh, colorscale, output_dir)
  % [lh_annot_file, rh_annot_file] = CBIG_ProjectFsaverageDataToDiscretizedAnnotation(lh_data, rh_data, discretization_res, min_thresh, max_thresh, output_dir)
  %
  % Project surface data to parcels on a surface annotation with the original values discretized
  %
  % Input:
  %   - lh_data           : a vector containing the surface data of the left hemisphere in either 'fsaverage5', 'fsaverage6', or 'fsaverage' template
  %   - rh_data           : a vector containing the surface data of the right hemisphere in either 'fsaverage5', 'fsaverage6', or 'fsaverage' template
  %   - discretization_res: number of distinct discrete values
  %   - min_thresh        : minimum threshold of the output values
  %   - max_thresh        : maximum threshold of the output values
  %   - colorscale        : a discrete colorscale in Matlab format
  %   - output_dir        : output directory
  % Output:
  %   - lh_annot_file     : containing the discretized surface annotation of the left hemisphere
  %   - rh_annot_file     : containing the discretized surface annotation of the right hemisphere
  %
  % Example:
  %   [lh_annot_file, rh_annot_file] = CBIG_ProjectFsaverageDataToDiscretizedAnnotation(lh_data, rh_data, 28, 1e-05, 5e-05, '/data/Work/area8_visualization')
  %   Discretize the surface data stored in lh_data and rh_data into 28 value levels between 1e-05 and 5e-05. The discretized datas are stored inside two annotation files lh_annot_file, rh_annot_file stored under '/data/Work/area8_visualization'
  %   
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  if length(lh_data) ~= length(rh_data)
    disp('Left and right hemispheres do not have the same number of data points!');
    exit;
  end
  
  if length(lh_data) == 40962 % fsaverage6
    ref_annot_dir = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage6', 'label');
  elseif length(lh_data) == 10242 % fsaverage5
    ref_annot_dir = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage5', 'label');
  elseif length(lh_data) == 163842 % fsaverage
    ref_annot_dir = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage', 'label');
  end
  
  [lh_vertices, lh_label, lh_colortable] = CBIG_MapSingleHemiSurfaceDataToAnnotation(lh_data, 'lh', ref_annot_dir, discretization_res, colorscale, min_thresh, max_thresh);
  lh_annot_file = fullfile(output_dir, 'lh_parcels_from_surf_data.annot');
  write_annotation(lh_annot_file, lh_vertices, lh_label, lh_colortable);
  [rh_vertices, rh_label, rh_colortable] = CBIG_MapSingleHemiSurfaceDataToAnnotation(rh_data, 'rh', ref_annot_dir, discretization_res, colorscale, min_thresh, max_thresh);
  rh_annot_file = fullfile(output_dir, 'rh_parcels_from_surf_data.annot');
  write_annotation(rh_annot_file, rh_vertices, rh_label, rh_colortable);
