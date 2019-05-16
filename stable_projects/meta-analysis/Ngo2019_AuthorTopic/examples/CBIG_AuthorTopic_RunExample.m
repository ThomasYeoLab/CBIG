function CBIG_AuthorTopic_RunExample
  % CBIG_AuthorTopic_RunExample
  %
  % Example of performing coordinate-based meta-analysis with the author-topic model
  % to discover cognitive components of self-generated thought
  % The inference was performed with the following hyperparameters
  %  - K = 2:3 - number of cognitive components to be estimated is 2 and 3.
  %    More components should be tried in actual experiments.
  %  - alpha = 100, eta = 0.01 - hyperparameters of the Dirichlet's distribution
  %  - seeds = 1:3 - 3 reinitializations per K.
  %    At least 1000 reinitializations should be used in actual experiments to
  %    achieve stable estimates.
  %  - Input data is activation foci of self-generated thought dataset saved at
  %    ../SelfGeneratedThought/MNI152_ActivationCoordinates/SelfGeneratedThought_AllCoordinates.txt
  %
  %  Output:
  %  - MNI152 activation coordinates are converted to CVB input format and
  %    saved at ./data/SelfGeneratedThought_CVBData.mat
  %  - Estimates of the author-topic model's parameters are saved at
  %    ./outputs/K<K>/alpha100_eta0.01/params_K<K>_SEED<seed>.mat
  %  - Visualization of the 2-component solution is saved under ./figures directory
  %  - BIC files and figures are saved under ./BIC directory
  %
  % Example:
  %   CBIG_AuthorTopic_RunExample
  %   Perfoming a coordinate-based meta-analysis with the author-topic model using self-generated thought data
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  tic
  
  % add paths to functions specific to author-topic model
  CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
  utilitiesDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
    'Ngo2019_AuthorTopic', 'utilities');
  addpath(fullfile(utilitiesDir, 'preprocessing'));
  addpath(fullfile(utilitiesDir, 'inference'));
  addpath(fullfile(utilitiesDir, 'BIC'));
  addpath(fullfile(utilitiesDir, 'visualization'));
  
  % prepare dataset
  textFilePath = fullfile(CBIG_CODE_DIR, 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', ...
    'SelfGeneratedThought', 'MNI152_ActivationCoordinates', 'SelfGeneratedThought_AllCoordinates.txt');
  dataDirPath = fullfile(pwd, 'data');
  dataFileName = 'SelfGeneratedThought_CVBData.mat';
  system(['mkdir -p ' dataDirPath]);
  CBIG_AuthorTopic_GenerateCVBDataFromText(textFilePath, dataDirPath, dataFileName);
  data = load(fullfile(dataDirPath, dataFileName));

  % inference
  allKs = 2:3; % set allKs to wider range number of components (e.g. 1:5) to discover the full range of solutions
  seeds = 1:2; % set seeds = 1:1000 to obtain stable estimates
  alpha = 100;
  eta = 0.01;
  doSmoothing = 1;
  workDir = pwd;
  cvbData = fullfile(dataDirPath, dataFileName);

  for K = allKs
    for seed = seeds
      CBIG_AuthorTopic_RunInference(seed, K, alpha, eta, doSmoothing, workDir, cvbData);
    end
  end

  % find best solutions from all reinitializations
  outputsDir = fullfile(workDir, 'outputs');
  for K = allKs
    CBIG_AuthorTopic_ComputeBestSolution(outputsDir, K, seeds, alpha, eta);
  end

  % visualize components, e.g. at 2-component solution
  figuresDir = fullfile(workDir, 'figures');
  bestSolutionDir = fullfile(workDir, 'outputs', 'bestSolution', ['alpha' num2str(alpha) '_eta' num2str(eta)]);
  inputFile = fullfile(bestSolutionDir, 'BestSolution_K002.mat');
  CBIG_AuthorTopic_VisualizeComponentsOnBrainSurface(inputFile, figuresDir);
  assert(exist(fullfile(figuresDir, 'clear_brain_min1e-5_max5e-5', 'C1.grid.png'), 'file') == 2, ...
    'Missing visualization for component C1');
  assert(exist(fullfile(figuresDir, 'clear_brain_min1e-5_max5e-5', 'C2.grid.png'), 'file') == 2, ...
    'Missing visualization for component C2');

  % compute BIC measures
  maskPath = fullfile(dataDirPath, 'mask', 'expMask.nii.gz');
  bicDir = fullfile(workDir, 'BIC');
  CBIG_AuthorTopic_ComputeBIC(allKs, maskPath, bestSolutionDir, bicDir);
  assert(exist(fullfile(bicDir, 'BIC_3dFWHMx.eps'), 'file') == 2, 'Missing BIC plot');

  % clean up paths
  rmpath(fullfile(utilitiesDir, 'preprocessing'));
  rmpath(fullfile(utilitiesDir, 'inference'));
  rmpath(fullfile(utilitiesDir, 'BIC'));
  rmpath(fullfile(utilitiesDir, 'visualization'));

  toc
