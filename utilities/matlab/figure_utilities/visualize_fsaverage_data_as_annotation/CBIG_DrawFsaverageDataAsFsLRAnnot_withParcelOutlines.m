function CBIG_DrawFsaverageDataAsFsLRAnnot_withParcelOutlines(lh_data, rh_data, abs_path_to_output_dir, label, colorscale, parcellation_name, discretization_res, min_thresh, max_thresh)

  % CBIG_DrawFsaverageDataAsFsLRAnnot_withParcelOutlines(lh_data, rh_data, abs_path_to_output_dir, label, colorscale, discretization_res, min_thresh, max_thresh)
  %
  % Wrapper function to draw fsaverage surface data stored in files as surface annotation and capture screenshots of their visualization in FreeView
  %
  % Input:
  %   Compulsory:
  %     - lh_data              : column vector containing the fsaverage surface data of the left hemisphere.
  %     - rh_data              : column vector containing the fsaverage surface data of the right hemisphere.
  %                              The vector can have length 10242 (fsaverage5), 40962 (fsaverage6), or 163842 (fsaverage)
  %     - abs_path_to_output_dir: absoluate path to the output folder.
  %     - label                : common label of the output files.
  %   Optional:
  %     - colorscale           : string specifying the colorscale.
  %                              'clear_brain' (default): Human Connectome Workbench's clear_brain color scheme
  %                              'parula': Matlab's parula color scheme
  %                              'hsv': Matlab's HSV color scheme
  %     - parcellation_name    : name of the reference parcellation, either 'aparc_annot' (default) or 'aparc_a2009s_annot'
  %     - discretization_res   : number of discrete intensity levels. Default: 28
  %     - min_thresh           : minimum threshold of the visualized data. Values below min_thresh are excluded. Default: 1e-5
  %     - max_thresh           : minimum threshold of the visualized data. Values below min_thresh are excluded. Default: 1
  %
  % Example:
  %   CBIG_DrawFsaverageDataAsAnnotationInFsLR(lh_data, rh_data, '/data/Work/area8_visualization', 'area8')
  %   Draw the surface data saved in lh_data and rh_data and save the resulting images under /data/Work/area8_visualization with the label area8 prefixed in all file names. The values in the resulting images are discretized in 28 levels, with the minimum at 1e-5 and maximum threshold at 1
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if nargin < 9
    max_thresh = 1;
  else
    if ischar(max_thresh)
      max_thresh = str2num(max_thresh);
    end
  end

  if nargin < 8
    min_thresh = 1e-5;
  else
    if ischar(min_thresh)
      min_thresh = str2num(min_thresh);
    end
  end

  if nargin < 7
    discretization_res = 28;
  else
    if ischar(discretization_res)
      discretization_res = str2num(discretization_res);
    end
  end
  
  if nargin < 6
    parcellation_name = 'aparc_annot';
  end
  
  if nargin < 5
    colorscale = CBIG_GenerateClearbrainColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
  elseif ~isnumeric(colorscale)
    if strcmp(colorscale, 'hsv')
      colorscale = CBIG_GenerateHSVColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
    elseif strcmp(colorscale, 'parula')
      colorscale = CBIG_GenerateParulaColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
    elseif strcmp(colorscale, 'clear_brain') || strcmp(colorscale, 'default')
      colorscale = CBIG_GenerateClearbrainColorscale(discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);
    end
  end

  tmp_dir = fullfile(abs_path_to_output_dir, 'tmp');
  system(['mkdir -p ' tmp_dir]);
  
  if size(lh_data, 1) == 1
    lh_data = lh_data';
  end
  
  if size(rh_data, 1) == 1
    rh_data = rh_data';
  end

  [lh_annot_file, rh_annot_file] = CBIG_ProjectFsaverageDataToDiscretizedAnnotation(lh_data, rh_data, discretization_res, min_thresh, max_thresh, colorscale, tmp_dir);

  parcel_dir = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'figure_utilities', 'visualize_fsaverage_data_as_annotation', 'parcellation_outlines');
  lh_parcel_outline_file = fullfile(parcel_dir, ['lh.' parcellation_name '.outline.annot']);
  rh_parcel_outline_file = fullfile(parcel_dir, ['rh.' parcellation_name '.outline.annot']);

  lh_combined_annot_file = fullfile(tmp_dir, 'lh.combined.annot');
  rh_combined_annot_file = fullfile(tmp_dir, 'rh.combined.annot');

  CBIG_CombineSurfaceAnnotations(lh_annot_file, lh_parcel_outline_file, lh_combined_annot_file);
  CBIG_CombineSurfaceAnnotations(rh_annot_file, rh_parcel_outline_file, rh_combined_annot_file);

  lh_final_annot_file = fullfile(tmp_dir, ['lh.' label '.annot']);
  rh_final_annot_file = fullfile(tmp_dir, ['rh.' label '.annot']);
  
  CBIG_AnnotateMedialWall(lh_combined_annot_file, 'lh', lh_final_annot_file);
  CBIG_AnnotateMedialWall(rh_combined_annot_file, 'rh', rh_final_annot_file);
 

  CBIG_VisualizeSurfaceAnnotationInFreeview(lh_final_annot_file, rh_final_annot_file, label, abs_path_to_output_dir);
