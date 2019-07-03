function CBIG_AuthorTopicEM_ReplicateBrainMapResults(seed, K, alpha, eta)
  % CBIG_AuthorTopicEM_ReplicateBrainMapResults
  %
  % Wrapper function to replicate the cognitive components estimated from BrainMap
  % from Yeo et al. 2015. Input data is stored in CBIG's HPC server.
  % 
  % The inference was performed with the following inputs
  %  - K: number of cognitive components to be estimated (to be set to 12)
  %  - alpha, eta: hyperparameters of the Dirichlet's distribution
  %    (default 100 and 0.01 respectively) 
  %  - seed: initialization seed of the random number's generator
  %
  %  Output:
  %  - Estimates of the author-topic model's parameters are saved at
  %    ./EM_outputs_GIBBS_init/K<K>/alpha100_eta0.01/seed<seed>.mat
  %
  % Example:
  %   CBIG_AuthorTopicEM_ReplicateBrainMapResults(42, 12, 100, 0.1)
  %   Perfominga coordinate-based meta-analysis with the author-topic model
  %   with K = 12 components, alpha = 100, eta = 0.01 and random number
  %   generator's initialized at 42.
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  % add paths to functions specific to author-topic model
  CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
  addpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
    'non_default_packages', 'topictoolbox')); % Gibbs sampler
  utilities_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
    'Yeo2015_AuthorTopicEM', 'utilities');
  addpath(fullfile(utilities_dir, 'inference'));
 
  brain_mask2mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'meta-analysis', 'Yeo2015_AuthorTopicEM', 'utilities', ...
    'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz'));

  % prepare dataset
  work_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
    'Yeo2015_AuthorTopicEM', 'replication');
  data_dir_path = fullfile('/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects', ...
    'meta-analysis/Yeo2015_AuthorTopicEM/data');
  gibbs_data_file_name = 'gibbs_AT_bin_unweighted_min1_max1.mat';
  em_data_file_name = 'EM_AT_bin_wc_unweighted_collapsed_min1_max1.mat';
  system(['mkdir -p ' data_dir_path]);
  
  % gibbs sampler's parameters
  gibbs_burn_in = 100; % set tto 100 or more in actual experiments
  gibbs_inverval = gibbs_burn_in;
  num_gibbs_samples = 1; % set to 100 or more in actual experiments
  gibbs_input_path = fullfile(data_dir_path, gibbs_data_file_name);

  disp('Start Gibbs sampling');
  CBIG_AuthorTopicEM_RunGibbsSampler(seed, K, gibbs_burn_in, gibbs_inverval, num_gibbs_samples, ...
      gibbs_input_path, work_dir, brain_mask2mm, alpha, eta);
  disp('Finished Gibbs sampling');

  em_init_type = 'GIBBS';
  em_input_path = fullfile(data_dir_path, em_data_file_name);
  disp('Start EM inference');
  gibbs_init_dir = fullfile(work_dir, 'GIBBS_outputs', ...
    ['K' num2str(K)], ['alpha' num2str(alpha) '_eta' num2str(eta)]);
  gibbs_init_file = fullfile(gibbs_init_dir, ['seed' num2str(seed)], ...
    ['seed' num2str(seed) '.burn' num2str(gibbs_burn_in) '.int' ...
    num2str(gibbs_inverval) '.burnin.mat']);
  CBIG_AuthorTopicEM_RunEM(seed, K, em_input_path, brain_mask2mm, ...
    alpha, eta, em_init_type, gibbs_init_file);
  disp('Finished EM inference');

  % clean up
  rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', ...
    'non_default_packages', 'topictoolbox')); % Gibbs sampler
  rmpath(fullfile(utilities_dir, 'inference'));
