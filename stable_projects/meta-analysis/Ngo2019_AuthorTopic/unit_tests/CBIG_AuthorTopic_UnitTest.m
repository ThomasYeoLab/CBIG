function CBIG_AuthorTopic_UnitTest
  % CBIG_AuthorTopic_UnitTest
  %
  % Wrapper function to unit-test a coordinate-based meta-analysis with the author-topic model.
  % The inference was performed with the following hard-coded hyperparameters
  %  - K = 1:3 - number of cognitive components to be estimated is set between 1 and 3 (inclusive).
  %  - alpha = 100, eta = 0.01 - hyperparameters of the Dirichlet's distribution
  %  - seeds = 1:3 - number of reinitializations per K is set between 1 and 3 (inclusive)
  %  - Input data is activation foci of self-generated thought dataset saved in ../SelfGeneratedThought/MNI152_ActivationCoordinates/SelfGeneratedThought_AllCoordinates.txt
  %
  %  Output:
  %  - Formatted input data is saved at ./data/SelfGeneratedThought_CVBData.mat
  %  - Estimates of the author-topic model's parameters are saved at ./outputs/K<K>/alpha100_eta0.01/params_K3_SEED<seed>.mat
  %  - Visualization of the 2-component solution is saved under ./figures directory
  %  - BIC files and figures are saved under ./BIC directory
  %
  % Example:
  %   CBIG_AuthorTopic_UnitTest
  %   Perfoming unit-testing of a coordinate-based meta-analysis with the author-topic model using self-generated thought data
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  % add paths to functions specific to author-topic model
  CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
  addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'preprocessing'));
  addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'inference'));
  addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'BIC'));
  addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'visualization'));
  
  % prepare dataset
  textFilePath = fullfile(CBIG_CODE_DIR, 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'SelfGeneratedThought', 'MNI152_ActivationCoordinates', 'SelfGeneratedThought_AllCoordinates.txt');
  dataDirPath = fullfile(pwd, 'data');
  dataFileName = 'SelfGeneratedThought_CVBData.mat';
  system(['mkdir -p ' dataDirPath]);
  CBIG_AuthorTopic_GenerateCVBDataFromText(textFilePath, dataDirPath, dataFileName);
  data = load(fullfile(dataDirPath, dataFileName));
  assert(size(data.act, 1) == 179, 'Wrong number of experiments');
  assert(size(data.act, 2) == 284100, 'Wrong number of voxels in the brain');
  assert(size(data.taskByExp, 1) == 7, 'Wrong number of tasks');
  assert(size(data.taskByExp, 2) == size(data.act, 1), 'Mismatched number of experiments');
  

  % inference
  allKs = 1:3; % set allKs to higher number of components to discover the full range of solutions
  alpha = 100;
  eta = 0.01;
  doSmoothing = 1;
  workDir = pwd;
  cvbData = fullfile(dataDirPath, dataFileName);
  seeds = 1:3; % set seeds = 1:1000 to obtain stable estimates

  for K = allKs
    for seed = seeds
      CBIG_AuthorTopic_RunInference(seed, K, alpha, eta, doSmoothing, workDir, cvbData);
    end
  end

  % find best solutions from all reinitializations
  outputsDir = fullfile(workDir, 'outputs');
  for K = allKs
    CBIG_AuthorTopic_ComputeBestSolution(outputsDir, K, seeds, alpha, eta);

    bestSolution = load(fullfile(workDir, 'outputs', 'bestSolution', ['alpha' num2str(alpha) '_eta' num2str(eta)], ['BestSolution_K' num2str(K, '%03d') '.mat']));
    assert(size(bestSolution.params.beta, 1) == K, 'Wrong number of componets in beta');
    assert(size(bestSolution.params.beta, 2) == size(data.act, 2), 'Mismatched number of voxels in beta');
    assert(size(bestSolution.params.theta, 2) == K, 'Wrong number of componets in theta');
    assert(size(bestSolution.params.theta, 1) == size(data.taskByExp, 1), 'Mismatched number of tasks in theta');
  end

  % visualize components, e.g. at 2-component solution
  figuresDir = fullfile(workDir, 'figures');
  bestSolutionDir = fullfile(workDir, 'outputs', 'bestSolution', ['alpha' num2str(alpha) '_eta' num2str(eta)]);
  inputFile = fullfile(bestSolutionDir, 'BestSolution_K002.mat');
  CBIG_AuthorTopic_VisualizeComponentsOnBrainSurface(inputFile, figuresDir);
  assert(exist(fullfile(figuresDir, 'clear_brain_min1e-5_max5e-5', 'C1.grid.png'), 'file') == 2, 'Missing visualization for component C1');
  assert(exist(fullfile(figuresDir, 'clear_brain_min1e-5_max5e-5', 'C2.grid.png'), 'file') == 2, 'Missing visualization for component C2');

  % compute BIC measures
  maskPath = fullfile(dataDirPath, 'mask', 'expMask.nii.gz');
  bicDir = fullfile(workDir, 'BIC');
  CBIG_AuthorTopic_ComputeBIC(allKs, maskPath, bestSolutionDir, bicDir);
  assert(exist(fullfile(bicDir, 'BIC_3dFWHMx.eps'), 'file') == 2, 'Missing BIC plot');

  % clean up
  rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'preprocessing'));
  rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'inference'));
  rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'BIC'));
  rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'visualization'));
