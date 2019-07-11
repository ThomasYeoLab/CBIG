function CBIG_AuthorTopic_VisualizeComponentsOnBrainSurface(paramsPath, ...
    outputDir, discretizationRes, minThresh, maxThresh, ...
    colorscaleName, doRemoveSmallClusters)
% CBIG_AuthorTopic_VisualizeComponentsOnBrainSurface(paramsPath, ...
%   outputDir, discretizationRes, minThresh, maxThresh, ...
%   colorscaleName, doCleanup)
%
% Visualize on the brain surface the probabilities of brain voxles being
% activated by the components (Pr(voxel | component)) of the author-topic
% model estimated by the Collapsed Variational Bayes algorithm.
%
% Input:
%  - paramsPath: absolute path to the .mat file containing the model
%                parameter estimates
%  - outputDir: absolute path to the directory containing the output
%               images and intermediate files. The final images are saved
%               under <outputDir>/<colorscaleName>_min<minThresh>_max<maxThresh>
%  - discretizationRes: the probabilities of brain voxles being activated
%               by components (Pr(voxel | component)) are discretized for
%               visualization. 'discretizationRes' determines the number of
%               discrete values used. Default value: 28
%  - minThresh: lower threshold of the values being visualized.
%               Default: 1e-5
%  - maxThresh: upper threshold of the values being visualized.
%               Default: 5e-5
%  - colorscaleName: name of the colorscale being used for visualization.
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
%   CBIG_AuthorTopic_VisualizeComponentsOnBrainSurface(...
%     '/Work/outputs/bestSolution/BestSolution_K002.mat', ...
%     '/Work/figures', 28, 1e-5, 5e-5, 'clear_brain', true);
%   Visualize the author-topic model's parameter estimate saved in
%   '/Work/outputs/bestSolution/BestSolution_K002.mat'. The final images
%   are saved under /Work/outputs/figures/clear_brain_min1e-5_max5e-5
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    load(paramsPath);
    brainMask = MRIread(params.maskPath);
    tmp = zeros([numel(brainMask.vol), params.K]);
    for i = 1:params.K
      tmp(brainMask.vol(:) == 1, i) = params.beta(i, :);
    end

    processes = brainMask;
    processes.nframes = params.K;
    processes.vol = reshape(tmp, [size(brainMask.vol) params.K]);

    if nargin < 7
      doRemoveSmallClusters = true;
    end

    if nargin < 6
      colorscaleName = 'clear_brain';
    end

    if nargin < 4
      minThresh = 1e-05;
      maxThresh = 5e-05;
      finalImagesDir = fullfile(outputDir, [colorscaleName '_min1e-5_max5e-5']);
    else
      finalImagesDir = fullfile(outputDir, [colorscaleName '_min' num2str(minThresh) '_max' num2str(maxThresh)]);
    end

    if nargin < 3
      discretizationRes = 28;
    end

    system(['mkdir -p ' finalImagesDir]);

    if strcmp(colorscaleName, 'hsv')
      colorscale = CBIG_GenerateHSVColorscale(discretizationRes, minThresh, maxThresh, finalImagesDir);
    elseif strcmp(colorscaleName, 'parula')
      colorscale = CBIG_GenerateParulaColorscale(discretizationRes, minThresh, maxThresh, finalImagesDir);
    elseif strcmp(colorscaleName, 'clear_brain') || strcmp(colorscaleName, 'default')
      colorscale = CBIG_GenerateClearbrainColorscale(discretizationRes, minThresh, maxThresh, finalImagesDir);
    end

    num_K = params.K;

    disp('Projecting to FSLR as ');
    [lhProjectedData, rhProjectedData] = CBIG_ProjectMNI2fsaverage_Ants(processes, 'fsaverage6');

    CBIG_AuthorTopic_VisualizeComponentsInFS_LR(num_K, finalImagesDir, ...
      discretizationRes, lhProjectedData, rhProjectedData, outputDir, minThresh, maxThresh, ...
      colorscale, doRemoveSmallClusters);

function CBIG_AuthorTopic_VisualizeComponentsInFS_LR(numK, finalImagesDir, ...
    discretization_res, lhProjectedData, rhProjectedData, fslrDir, minThresh, maxThresh, ...
    colorscale, doRemoveSmallClusters)
  SMOOTH = 'metric';

  SURF_COMP_SIZE_THRESH = 20;

  lhAvgMesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k', 'inflated', 'aparc.annot');
  rhAvgMesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k', 'inflated', 'aparc.annot');

  [lhVertices, refLhLabels, refLhColortable] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), ...
  'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'lh.aparc.annot'));
  [rhVertices, refRhLabels, refRhColortable] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), ...
  'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'rh.aparc.annot'));

  for K = 1:numK
    componentDir = fullfile(fslrDir, ['C' num2str(K)]);


    if exist(componentDir, 'dir')
      system(['rm -r ' componentDir]);
    end

    lhProjectedComponent = lhProjectedData(K, :)';
    rhProjectedComponent = rhProjectedData(K, :)';
    disp('Transformation with wb_command');

    [origLhFslr32kProjectedComponent, origRhFslr32kProjectedComponent, ~, ~] = CBIG_project_fsaverage2fsLR(...
      lhProjectedComponent, rhProjectedComponent, 'fsaverage6', SMOOTH, componentDir);
    [origLhLabels, lhColortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(...
      origLhFslr32kProjectedComponent, discretization_res, colorscale, minThresh, maxThresh);
    [origRhLabels, rhColortable] = CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(...
      origRhFslr32kProjectedComponent, discretization_res, colorscale, minThresh, maxThresh);

    num_colors = size(lhColortable.table, 1);
    underlayLabels = lhColortable.table(num_colors, 5);
    refMedialwallLabels = refLhColortable.table(1, 5);

    binaryLhLabels = zeros(size(origLhLabels));
    binaryRhLabels = zeros(size(origRhLabels));
    binaryLhLabels(origLhLabels ~= underlayLabels) = 1;
    binaryRhLabels(origRhLabels ~= underlayLabels) = 1;

    if (doRemoveSmallClusters)
      lhLabelMask = CBIG_RemoveIsolatedSurfaceComponents(lhAvgMesh, binaryLhLabels, SURF_COMP_SIZE_THRESH);
      rhLabelMask = CBIG_RemoveIsolatedSurfaceComponents(rhAvgMesh, binaryRhLabels, SURF_COMP_SIZE_THRESH);
    end

    lhFslr32kLabels = origLhLabels;
    lhFslr32kLabels(lhLabelMask == 0) = underlayLabels;

    rhFslr32kLabels = origRhLabels;
    rhFslr32kLabels(rhLabelMask == 0) = underlayLabels;

    mkdir(componentDir);
    lhFslr32kAnnotFile = fullfile(componentDir, 'lh_fslr_parcels.annot');
    rhFslr32kAnnotFile = fullfile(componentDir, 'rh_fslr_parcels.annot');

    CBIG_AnnotateSingleHemiMedialWall(lhVertices, lhFslr32kLabels, lhColortable, ...
      refLhLabels, refMedialwallLabels, lhFslr32kAnnotFile);
    CBIG_AnnotateSingleHemiMedialWall(rhVertices, rhFslr32kLabels, rhColortable, ...
      refRhLabels, refMedialwallLabels, rhFslr32kAnnotFile);

    CBIG_VisualizeSurfaceAnnotationInFreeview(lhFslr32kAnnotFile, rhFslr32kAnnotFile, ...
      'fs_LR_32k', ['C' num2str(K)], finalImagesDir);
  end
  close all;
