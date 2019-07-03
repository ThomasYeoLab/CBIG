function CBIG_AuthorTopicEM_ReplicateMostTypicalOrBestSolutionWrapper
  % CBIG_AuthorTopicEM_RunEMwithGibbsInitExample
  %
  % Wrapper function to replicate the most typical and best1 12-component solution
  %  Output:
  %  - Visualization of the most typical/best solution is saved under ./figures directory
  %
  % Example:
  %   CBIG_AuthorTopicEM_ReplicateMostTypicalOrBestSolutionWrapper
  %   Compute the most typical and best 12-component solution
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  K = 12;
  SEEDS = 1:100;
  ALPHA = 100;
  ETA = 0.01;

  work_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
    'Yeo2015_AuthorTopicEM', 'replication');
  em_init_type = 'GIBBS';
  output_dir = fullfile(work_dir, ['EM_outputs_' em_init_type '_init']);
  
  % add paths to functions specific to author-topic model
  CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
  utilities_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
    'Yeo2015_AuthorTopicEM', 'utilities');
  addpath(fullfile(utilities_dir, 'inference'));
  addpath(fullfile(utilities_dir, 'visualization'));

  % find the best solution from all reinitializations
  disp('Compute the best solution');
  best_solution_dir = fullfile(output_dir, 'best_solution', ['alpha' num2str(ALPHA) '_eta' num2str(ETA)]);
  CBIG_AuthorTopicEM_ComputeBestSolution(output_dir, K, SEEDS, ALPHA, ETA);

  % find the most typical solutions from all reinitializations
  disp('Compute the most typical solution');
  output_dir = fullfile(work_dir, ['EM_outputs_' em_init_type '_init']);
  avg_solution_dir = fullfile(output_dir, 'avg_solution', ['alpha' num2str(ALPHA) '_eta' num2str(ETA)]);
  CBIG_AuthorTopicEM_ComputeMostTypicalSolution(output_dir, K, SEEDS, ALPHA, ETA);

  % visualize most typical solution
  figures_dir = fullfile(work_dir, 'figures');
  input_path = fullfile(avg_solution_dir, ['avg_solution_K' num2str(K, '%03d') '.mat']);
  CBIG_AuthorTopicEM_VisualizeComponentsOnBrainSurface(input_path, figures_dir);

  % clean up
  rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
    'non_default_packages', 'topictoolbox')); % Gibbs sampler
  rmpath(fullfile(utilities_dir, 'inference'));
  rmpath(fullfile(utilities_dir, 'visualization'));
