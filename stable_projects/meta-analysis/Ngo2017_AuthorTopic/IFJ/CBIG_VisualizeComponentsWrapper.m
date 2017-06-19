function CBIG_VisualizeComponentsWrapper(K)
  ALPHA = 100;
  ETA = 0.01;

  currDir = pwd;
  alphaEtaPrefix = ['alpha' num2str(ALPHA) '_eta' num2str(ETA)];

  solutionPath = fullfile(currDir, 'outputs', 'best_solution', alphaEtaPrefix, ['BestSolution_K' num2str(K, '%03d') '.mat']);
  outputDir = fullfile(currDir, 'outputs', 'best_solution', alphaEtaPrefix, 'figures');
  minThresh = 1e-05;
  maxThresh = 5e-05;
  colorscale = 'clear_brain';

  CBIG_PlotComponentsToSurfaceAsParcels_withCleanUp(solutionPath, outputDir, minThresh, maxThresh, colorscale);

  cd(currDir);
