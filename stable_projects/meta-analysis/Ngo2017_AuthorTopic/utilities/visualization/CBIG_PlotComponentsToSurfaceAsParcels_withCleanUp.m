function CBIG_PlotComponentsToSurfaceAsParcels_withCleanUp(input_path, fslr_dir, min_thresh, max_thresh, colorscale)
    % input_path - the input mat file. The brain activation is inside params.beta
    % flsr_dir - output for FS_LR template directory

    load(input_path);
    brain_mask = MRIread('~ngohgia/templates/MNI_mask_conformed.2mm.0.1.nii.gz');
    tmp = zeros([size(brain_mask.vol, 1)*size(brain_mask.vol, 2)*size(brain_mask.vol, 3), params.K]);
    for i = 1:params.K
      tmp(brain_mask.vol(:) == 1, i) = params.beta(i, :);
    end

    processes = brain_mask;
    processes.nframes = params.K;
    processes.vol = reshape(tmp, [size(brain_mask.vol) params.K]);
    
    if nargin < 5
      colorscale = 'clear_brain';
    end
    
    if nargin < 3
      min_thresh = ones(params.K, 1) * 1e-05;
      max_thresh = ones(params.K, 1) * 5e-05;
      final_images_dir = fullfile(fslr_dir, [colorscale '_min' num2str(min_thresh) '_max' num2str(max_thresh) '_clean']);     
    else
      if ~isnumeric(max_thresh)
        final_images_dir = fullfile(fslr_dir, [colorscale '_min' num2str(1e-05) '_maxVariable_clean']);          
        min_thresh = ones(params.K, 1) * 1e-05;
        max_thresh = max(params.beta, [], 2);    
      elseif length(max_thresh) == 1
        final_images_dir = fullfile(fslr_dir, [colorscale '_min' num2str(min_thresh) '_max' num2str(max_thresh) '_clean']);
        min_thresh = ones(params.K, 1) * min_thresh;
        max_thresh = ones(params.K, 1) * max_thresh;
      end
    end
    
    system(['mkdir -p ' final_images_dir]);

    num_K = params.K;
    
    disp('Projecting to FSLR as ');
    [lh_proj_data, rh_proj_data] = CBIG_ProjectMNI2fsaverage(processes, 'fsaverage6');
    
    ProjectActivationAsParcelsToFS_LR(num_K, final_images_dir, lh_proj_data, rh_proj_data, fslr_dir, min_thresh, max_thresh, colorscale);

function ProjectActivationAsParcelsToFS_LR(num_K, final_images_dir, lh_proj_data, rh_proj_data, fslr_dir, min_thresh, max_thresh, colorscale)
  SMOOTH = 'METRIC_AVERAGE_TILE';
  DISCRETIZATION_RES = 28;
  
  SURF_COMP_SIZE_THRESH = 20;
    
  lh_avg_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k', 'inflated', 'aparc.annot');
  rh_avg_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k', 'inflated', 'aparc.annot');
      
  for K = 1:num_K
    component_dir = fullfile(fslr_dir, ['C' num2str(K) '_clean']);
    if ~exist(component_dir, 'dir')
      mkdir(component_dir);
    end
      
    component_lh_proj = lh_proj_data(K, :)';
    component_rh_proj = rh_proj_data(K, :)';
    disp('Transformation with wb_command');
    
    [orig_lh_proj_32k, orig_proj_32k, ~, ~] = CBIG_project_from_fsaverage_to_fslr(component_lh_proj, component_rh_proj, component_dir, SMOOTH);
    [tmp_lh_annot_file, tmp_rh_annot_file] = CBIG_DiscretizeSurfDataAsParcels(orig_lh_proj_32k, orig_proj_32k, DISCRETIZATION_RES, min_thresh(K), max_thresh(K), component_dir, colorscale);
    [lh_vertices, orig_lh_labels, lh_color_table] = read_annotation(tmp_lh_annot_file);
    [rh_vertices, orig_rh_labels, rh_color_table] = read_annotation(tmp_rh_annot_file);
    
    num_colors = size(lh_color_table.table, 1);
    underlay_label = lh_color_table.table(num_colors-1, 5);
    medial_wall_label = lh_color_table.table(num_colors, 5);
    
    bin_lh_labels = zeros(size(orig_lh_labels));
    bin_rh_labels = zeros(size(orig_rh_labels));
    bin_lh_labels(orig_lh_labels ~= underlay_label) = 1;
    bin_rh_labels(orig_rh_labels ~= underlay_label) = 1;
    
    lh_label_mask = CBIG_RemoveIsolatedSurfaceComponents(lh_avg_mesh, bin_lh_labels, SURF_COMP_SIZE_THRESH);
    rh_label_mask = CBIG_RemoveIsolatedSurfaceComponents(rh_avg_mesh, bin_rh_labels, SURF_COMP_SIZE_THRESH);
    
    lh_32k_labels = orig_lh_labels;
    lh_32k_labels(lh_label_mask == 0) = underlay_label;
    lh_32k_labels(orig_lh_labels == medial_wall_label) = medial_wall_label;
    rh_32k_labels = orig_rh_labels;
    rh_32k_labels(rh_label_mask == 0) = underlay_label;
    rh_32k_labels(orig_rh_labels == medial_wall_label) = medial_wall_label;
    
    lh_fslr_annot_file = fullfile(component_dir, 'lh_fslr_parcels.annot');
    write_annotation(lh_fslr_annot_file, lh_vertices, lh_32k_labels, lh_color_table)
    rh_fslr_annot_file = fullfile(component_dir, 'rh_fslr_parcels.annot');
    write_annotation(rh_fslr_annot_file, rh_vertices, rh_32k_labels, rh_color_table)
    
    
    PlotAnnotationInFreeview(K, component_dir, final_images_dir, lh_fslr_annot_file, rh_fslr_annot_file);
    
    % system(['rm ' fullfile(component_dir, '*.tmp.*')]);
  end
  % close all;
  
function PlotAnnotationInFreeview(K, component_dir, final_images_dir, lh_annot_file, rh_annot_file)
  cd(component_dir);
  system('mkdir -p fv_output');
  
  colorscale_file = fullfile(component_dir, 'colorscale.png');
  system(CBIG_TrimFreeviewScreenshotCommand(colorscale_file, colorscale_file, 0.1, 0.40));
  
  raw_lh_lateral_file = fullfile('fv_output', ['C' num2str(K) '_lh_lateral.raw.png']);
  lh_lateral_file = fullfile('fv_output', ['C' num2str(K) '_lh_lateral.jpg']);
  tmp_lh_lateral_file = fullfile('fv_output', ['C' num2str(K) '_lh_lateral.tmp.png']);
  system(ProjectLHAnnotation(lh_annot_file, raw_lh_lateral_file, 0));
  system(ReplaceBlackByWhiteBkgCommand(raw_lh_lateral_file, tmp_lh_lateral_file));
  
  raw_lh_medial_file = fullfile('fv_output', ['C' num2str(K) '_lh_medial.raw.png']);
  lh_medial_file = fullfile('fv_output', ['C' num2str(K) '_lh_medial.jpg']);
  tmp_lh_medial_file = fullfile('fv_output', ['C' num2str(K) '_lh_medial.tmp.png']);
  system(ProjectLHAnnotation(lh_annot_file, raw_lh_medial_file, 180));
  system(ReplaceBlackByWhiteBkgCommand(raw_lh_medial_file, tmp_lh_medial_file));
  
  raw_rh_lateral_file = fullfile('fv_output', ['C' num2str(K) '_rh_lateral.raw.png']);
  rh_lateral_file = fullfile('fv_output', ['C' num2str(K) '_rh_lateral.jpg']);
  tmp_rh_lateral_file = fullfile('fv_output', ['C' num2str(K) '_rh_lateral.tmp.png']);
  system(ProjectRHAnnotation(rh_annot_file, raw_rh_lateral_file, 180));
  system(ReplaceBlackByWhiteBkgCommand(raw_rh_lateral_file, tmp_rh_lateral_file));
  
  raw_rh_medial_file = fullfile('fv_output', ['C' num2str(K) '_rh_medial.raw.png']);
  rh_medial_file = fullfile('fv_output', ['C' num2str(K) '_rh_medial.jpg']);
  tmp_rh_medial_file = fullfile('fv_output', ['C' num2str(K) '_rh_medial.tmp.png']);
  system(ProjectRHAnnotation(rh_annot_file, raw_rh_medial_file, 0));
  system(ReplaceBlackByWhiteBkgCommand(raw_rh_medial_file, tmp_rh_medial_file));
  
  tmp_grid_output_path = ['C' num2str(K) '_grid.raw.png'];
  grid_output_path = fullfile(final_images_dir, ['C' num2str(K) '_grid.png']);
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_lateral_file, lh_lateral_file));
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_medial_file, lh_medial_file));
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_lateral_file, rh_lateral_file));
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_medial_file, rh_medial_file));
  system(CBIG_CreateImagesGridCommand(tmp_grid_output_path, lh_lateral_file, lh_medial_file, rh_lateral_file, rh_medial_file));
  system(CBIG_InsertColorscaleToImage(colorscale_file, tmp_grid_output_path, grid_output_path));
  % system(CompressImgCommand(tmp_grid_output_path, grid_output_path));

  tmp_series_output_path = ['C' num2str(K) '_series.raw.png'];
  series_output_path = fullfile(final_images_dir, ['C' num2str(K) '_series.png']);
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_lateral_file, lh_lateral_file, 0.3, 0));
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_lh_medial_file, lh_medial_file, 0.3, 0));
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_lateral_file, rh_lateral_file, 0.3, 0));
  system(CBIG_TrimFreeviewScreenshotCommand(tmp_rh_medial_file, rh_medial_file, 0.3, 0));
  system(CBIG_CreateImagesSeriesCommand(tmp_series_output_path, lh_lateral_file, lh_medial_file, rh_lateral_file, rh_medial_file));
  system(CBIG_InsertColorscaleToImage(colorscale_file, tmp_series_output_path, series_output_path, 0.6, 0.82, 0.85));
  % system(CompressImgCommand(tmp_series_output_path, series_output_path));
  
  cd(pwd);

function command = ProjectLHAnnotation(annotation_file, output_file, degree)
  underlay = fullfile(getenv('CBIG_CODE_DIR'), 'data/templates/surface/fs_LR_32k/surf/lh.inflated');
  command = ['freeview -f ' underlay ':annotation=' annotation_file ...
      ':edgethickness=0:color=lightgray' ...
      ' -viewport 3d' ...
      ' -cam azimuth ' num2str(degree) ...
      ' -ss ' output_file ];
  
function command = ProjectRHAnnotation(annotation_file, output_file, degree)
  underlay = fullfile(getenv('CBIG_CODE_DIR'), 'data/templates/surface/fs_LR_32k/surf/rh.inflated');
  command = ['freeview -f ' underlay ':annotation=' annotation_file ...
      ':edgethickness=0:color=lightgray' ...
      ' -viewport 3d' ...
      ' -cam azimuth ' num2str(degree) ...
      ' -ss ' output_file ];
  
function command = CBIG_InsertColorscaleToImage(colorscale_file, underlay_image, final_image, new_dim_pct, x_offset_pct, y_offset_pct)
  if nargin < 4  
    new_dim_pct = 0.5;
    x_offset_pct = 0.33;
    y_offset_pct = 0.47;
  end
  
  if isnumeric(new_dim_pct)
    new_dim_pct = num2str(new_dim_pct);
  end
  
  if isnumeric(x_offset_pct)
    x_offset_pct = num2str(x_offset_pct);
  end
  
  if isnumeric(y_offset_pct)
    y_offset_pct = num2str(y_offset_pct);
  end
  
  small_colorscale_file = [colorscale_file '.resized.png'];
  
  command = ['w=`identify ' colorscale_file '| cut -f 3 -d " " | sed s/x.*//`; ' ...
      'new_w=`awk "BEGIN {printf(\"%.0f\", $w * ' new_dim_pct ')}"`; ' ...
      'underlay_w=`identify ' underlay_image '| cut -f 3 -d " " | sed s/x.*//`; ' ...
      'colorscale_x=`awk "BEGIN {printf(\"%.0f\", $underlay_w * ' x_offset_pct ')}"`; ' ...
      'h=`identify ' colorscale_file '| cut -f 3 -d " " | sed s/.*x//`; ' ...
      'new_h=`awk "BEGIN {printf(\"%.0f\", $h * ' new_dim_pct ')}"`; ' ...
      'underlay_h=`identify ' underlay_image '| cut -f 3 -d " " | sed s/.*x//`; ' ...
      'colorscale_y=`awk "BEGIN {printf(\"%.0f\", $underlay_h * ' y_offset_pct ')}"`; ' ...
      'convert ' colorscale_file ' -resize ${new_w}x${new_h} ' small_colorscale_file '; ' ...
      'composite -compose atop -geometry +${colorscale_x}+${colorscale_y} ' small_colorscale_file ' ' underlay_image ' ' final_image];  
