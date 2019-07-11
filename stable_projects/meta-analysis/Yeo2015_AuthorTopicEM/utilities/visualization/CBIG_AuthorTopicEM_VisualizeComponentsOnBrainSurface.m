function CBIG_AuthorTopicEM_VisualizeComponentsOnBrainSurface(params_path, ...
    output_dir, discretization_res, min_thresh, max_thresh, ...
    colorscale_name, do_remove_small_clusters)
% CBIG_AuthorTopicEM_VisualizeComponentsOnBrainSurface(params_path, ...
%   output_dir, discretization_res, min_thresh, max_thresh, ...
%   colorscale_name, do_cleanup)
%
% Visualize on the brain surface the probabilities of brain voxles being
% activated by the components (Pr(voxel | component)) of the author-topic
% model estimated by the Expectation-Maximization algorithm.
%
% Input:
%  - params_path: absolute path to the .mat file containing the model
%                parameter estimates
%  - output_dir: absolute path to the directory containing the output
%               images and intermediate files. The final images are saved
%               under <output_dir>/<colorscale_name>_min<min_thresh>_max<max_thresh>
%  - discretization_res: the probabilities of brain voxles being activated
%               by components (Pr(voxel | component)) are discretized for
%               visualization. 'discretization_res' determines the number of
%               discrete values used. Default value: 28
%  - min_thresh: lower threshold of the values being visualized.
%               Default: 1e-5
%  - max_thresh: upper threshold of the values being visualized.
%               Default: 5e-5
%  - colorscale_name: name of the colorscale being used for visualization.
%    Possible values are:
%    + clear_brain: Human Connectome Workbench's clear_brain color palette.
%    + hsv: Matlab's HSV colorscale.
%    + parula: Matlab's parula colorscale.
%    Default: clear_brain
%  - doRemoveSmallClusters: if true, clusters of voxels of fewer than 20 voxels are
%    removed from the visualization to produce a clearer presentation of
%    the most dominant patterns. Default: true
%
% Example:
%   CBIG_AuthorTopicEM_VisualizeComponentsOnBrainSurface(...
%     '/Work/outputs/bestSolution/BestSolution_K002.mat', ...
%     '/Work/figures', 28, 1e-5, 5e-5, 'clear_brain', true);
%   Visualize the author-topic model's parameter estimate saved in
%   '/Work/outputs/bestSolution/BestSolution_K002.mat'. The final images
%   are saved under /Work/outputs/figures/clear_brain_min1e-5_max5e-5
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    load(params_path);
    tmp = zeros([numel(params.brain_mask.vol), params.K]);
    for i = 1:params.K
      tmp(params.brain_mask.vol(:) == 1, i) = params.beta(i, :);
    end

    processes = params.brain_mask;
    processes.nframes = params.K;
    processes.vol = reshape(tmp, [size(params.brain_mask.vol) params.K]);

    if nargin < 7
      do_remove_small_clusters = true;
    end

    if nargin < 6
      colorscale_name = 'clear_brain';
    end

    if nargin < 4
      min_thresh = 1e-05;
      max_thresh = 5e-05;
      final_images_dir = fullfile(output_dir, [colorscale_name '_min1e-5_max5e-5']);
    else
      final_images_dir = fullfile(output_dir, [colorscale_name '_min' num2str(min_thresh) '_max' num2str(max_thresh)]);
    end

    if nargin < 3
      discretization_res = 28;
    end

    system(['mkdir -p ' final_images_dir]);

    if strcmp(colorscale_name, 'hsv')
      colorscale = CBIG_GenerateHSVColorscale(discretization_res, min_thresh, max_thresh, final_images_dir);
    elseif strcmp(colorscale_name, 'parula')
      colorscale = CBIG_GenerateParulaColorscale(discretization_res, min_thresh, max_thresh, final_images_dir);
    elseif strcmp(colorscale_name, 'clear_brain') || strcmp(colorscale_name, 'default')
      colorscale = CBIG_GenerateClearbrainColorscale(discretization_res, min_thresh, max_thresh, final_images_dir);
    end

    num_K = params.K;

    disp('Projecting to FSLR as ');
    [lh_projected_data, rh_projected_data] = CBIG_ProjectMNI2fsaverage_Ants(processes, 'fsaverage6');

    CBIG_AuthorTopic_VisualizeComponentsInFS_LR(num_K, final_images_dir, ...
      discretization_res, lh_projected_data, rh_projected_data, output_dir, min_thresh, max_thresh, ...
      colorscale, do_remove_small_clusters);

function CBIG_AuthorTopic_VisualizeComponentsInFS_LR(numK, final_images_dir, ...
    discretization_res, lh_projected_data, rh_projected_data, fslr_dir, min_thresh, max_thresh, ...
    colorscale, do_remove_small_clusters)
  SMOOTH = 'metric';

  SURF_COMP_SIZE_THRESH = 20;

  lh_avg_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k', 'inflated', 'aparc.annot');
  rh_avg_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k', 'inflated', 'aparc.annot');

  [lh_vertices, ref_lh_labels, ref_lh_colortable] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), ...
    'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'lh.aparc.annot'));
  [rh_vertices, ref_rh_labels, ref_rh_colortable] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), ...
    'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'rh.aparc.annot'));

  for K = 1:numK
    component_dir = fullfile(fslr_dir, ['C' num2str(K)]);

    if exist(component_dir, 'dir')
      system(['rm -r ', component_dir]);
    end

    lh_projected_component = lh_projected_data(K, :)';
    rh_projected_component = rh_projected_data(K, :)';
    disp('Transformation with wb_command');

    [orig_lh_fslr_32k_projected_component, orig_rh_fslr32k_projected_component, ~, ~] = ...
      CBIG_project_fsaverage2fsLR(lh_projected_component, rh_projected_component, 'fsaverage6', SMOOTH, component_dir);
  
    [orig_lh_labels, lh_colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation( ...
      orig_lh_fslr_32k_projected_component, discretization_res, colorscale, min_thresh, max_thresh);
    [orig_rh_labels, rh_colortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation( ...
      orig_rh_fslr32k_projected_component, discretization_res, colorscale, min_thresh, max_thresh);

    num_colors = size(lh_colortable.table, 1);
    underlay_labels = lh_colortable.table(num_colors, 5);
    ref_medialwall_labels = ref_lh_colortable.table(1, 5);

    binary_lh_labels = zeros(size(orig_lh_labels));
    binary_rh_labels = zeros(size(orig_rh_labels));
    binary_lh_labels(orig_lh_labels ~= underlay_labels) = 1;
    binary_rh_labels(orig_rh_labels ~= underlay_labels) = 1;

    if (do_remove_small_clusters)
      lh_label_mask = CBIG_RemoveIsolatedSurfaceComponents(lh_avg_mesh, binary_lh_labels, SURF_COMP_SIZE_THRESH);
      rh_label_mask = CBIG_RemoveIsolatedSurfaceComponents(rh_avg_mesh, binary_rh_labels, SURF_COMP_SIZE_THRESH);
    end

    lh_fslr32k_labels = orig_lh_labels;
    lh_fslr32k_labels(lh_label_mask == 0) = underlay_labels;

    rh_fslr32k_labels = orig_rh_labels;
    rh_fslr32k_labels(rh_label_mask == 0) = underlay_labels;

    mkdir(component_dir)
    lh_fslr32k_annotFile = fullfile(component_dir, 'lh_fslr_parcels.annot');
    rh_fslr32k_annotFile = fullfile(component_dir, 'rh_fslr_parcels.annot');

    CBIG_AnnotateSingleHemiMedialWall(lh_vertices, lh_fslr32k_labels, lh_colortable, ...
      ref_lh_labels, ref_medialwall_labels, lh_fslr32k_annotFile);
    CBIG_AnnotateSingleHemiMedialWall(rh_vertices, rh_fslr32k_labels, rh_colortable, ...
      ref_rh_labels, ref_medialwall_labels, rh_fslr32k_annotFile);

    CBIG_VisualizeSurfaceAnnotationInFreeview(lh_fslr32k_annotFile, rh_fslr32k_annotFile, ...
      'fs_LR_32k', ['C' num2str(K)], final_images_dir);
  
    close all;
  end
