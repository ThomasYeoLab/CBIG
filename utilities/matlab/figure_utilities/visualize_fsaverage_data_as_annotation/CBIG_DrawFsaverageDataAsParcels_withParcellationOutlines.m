function CBIG_DrawFsaverageDataAsParcels_withParcellationOutlines(lh_data, rh_data, abs_path_to_output_dir, label, parcellation_name, discretization_res, min_thresh, max_thresh)
  % CBIG_DrawFsaverageDataAsParcels_withParcellationOutlines(lh_data, rh_data, abs_path_to_output_dir, label, min_thresh)
  %
  % Draw surface data in fsaverage template as parcels and visualize on the fs_LR with FreeView. Parcellation outlines are superimposed on the input fsaverage data for reference.
  %
  % Input:
  %   Compulsory:
  %     - lh_data              : left hemisphere surface data.
  %     - rh_data              : right hemisphere surface data.
  %                              The data can be in fsaverage5, fsaverage6 and fsaverage surface template.
  %     - ab_path_to_output_dir: absoluate path to the output folder.
  %     - label                : common label of the output files.
  %   Optional:
  %     - parcellation_name    : name of the reference parcellation, either 'aparc_annot' (default) or 'aparc_a2009s_annot'
  %     - discretization_res   : number of discrete intensity levels. Default: 28
  %     - min_thresh           : minimum threshold of the visualized data. Values below min_thresh are excluded. Default: 1e-5
  %     - max_thresh           : minimum threshold of the visualized data. Values below min_thresh are excluded. Default: 1
  %
  % Example:
  %   CBIG_DrawFsaverageDataAsParcels(lh_data, rh_data, '/data/Work/area8_visualization', 'area8', 1e-5)
  
  if nargin < 5
    parcellation_name = 'aparc_annot';
  end
 
  if nargin < 6
    discretization_res = 28;
  else
    if ischar(discretization_res)
      discretization_res = str2num(discretization_res);
    end
  end

  if nargin < 7
    min_thresh = 1e-5;
  else
    if ischar(min_thresh)
      min_thresh = str2num(min_thresh);
    end
  end

  if nargin < 8
    max_thresh = 1;
  else
    if ischar(max_thresh)
      max_thresh = str2num(max_thresh);
    end
  end

  [lh_annot_file, rh_annot_file] = CBIG_DiscretizeSurfaceDataAsParcels(lh_data', rh_data', discretization_res, min_thresh, max_thresh, abs_path_to_output_dir);

  lh_parcel_outline_file = fullfile(pwd, 'parcellation_outlines', ['lh.' parcellation_name '.outline.annot']);
  rh_parcel_outline_file = fullfile(pwd, 'parcellation_outlines', ['rh.' parcellation_name '.outline.annot']);

  lh_combined_annot_file = fullfile(abs_path_to_output_dir, 'lh.combined.annot');
  lh_final_annot_file = fullfile(abs_path_to_output_dir, ['lh.' label '.annot']);
  rh_combined_annot_file = fullfile(abs_path_to_output_dir, 'rh.combined.annot');
  rh_final_annot_file = fullfile(abs_path_to_output_dir, ['rh.' label '.annot']);
  
  CBIG_CombineSurfaceAnnotations(lh_annot_file, lh_parcel_outline_file, lh_combined_annot_file);
  CBIG_AnnotateMedialWall(lh_combined_annot_file, 'lh', lh_final_annot_file);
  
  CBIG_CombineSurfaceAnnotations(rh_annot_file, rh_parcel_outline_file, rh_combined_annot_file);
  CBIG_AnnotateMedialWall(rh_combined_annot_file, 'rh', rh_final_annot_file);
  
  CBIG_VisualizeSurfaceAnnotationInFreeview(abs_path_to_output_dir, lh_final_annot_file, rh_final_annot_file, label);
