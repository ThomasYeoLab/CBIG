classdef CBIG_AuthorTopicEM_unit_test < matlab.unittest.TestCase
  % CBIG_AuthorTopicEM_UnitTest
  %
  % Wrapper function to unit-test a coordinate-based meta-analysis with
  % the author-topic model and expectation-maximization (EM) algorithm.
  % The inference was performed with the following hard-coded hyperparameters
  %  - K = 2 - number of cognitive components to be estimated is set to 2
  %  - alpha = 100, eta = 0.01 - hyperparameters of the Dirichlet's distribution
  %  - seeds = 1:3 - number of reinitializations per K is set between 1 and 3 (inclusive)
  %  - Input data is activation foci of self-generated thought dataset saved in
  %    unit_test_sample_coordinates.tx
  %
  %  Output:
  %  - Formatted input data is saved at ./data/
  %  - Estimates of the author-topic model's parameters are saved at EM_outputs_GIBBS_init
  %  - Visualization of the 2-component solution is saved under ./figures directory
  %
  % Example:
  %   CBIG_AuthorTopicEM_UnitTest
  %   Perfoming unit-testing of a coordinate-based meta-analysis with the author-topic
  %   model using self-generated thought data
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  methods (Test)
    function test_authorTopic_EM(testCase)
      % add paths to functions specific to author-topic model
      CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
      addpath(fullfile(CBIG_CODE_DIR, 'external_packages', 'matlab', ...
        'non_default_packages', 'topictoolbox')); % Gibbs sampler
      utilities_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'meta-analysis', ...
        'Yeo2015_AuthorTopicEM', 'utilities');
      addpath(fullfile(utilities_dir, 'preprocessing'));
      addpath(fullfile(utilities_dir, 'inference'));
      addpath(fullfile(utilities_dir, 'visualization'));
     
      brain_mask1mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'data', ...
        'templates', 'volume','FSL_MNI152_FS4.5.0', 'mri', 'brain.mgz'));
      brain_mask2mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
        'meta-analysis', 'Yeo2015_AuthorTopicEM', 'utilities', ...
        'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz'));
    
      % prepare dataset
      work_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', ...
        'Yeo2015_AuthorTopicEM', 'unit_tests');
      text_data_path = fullfile(work_dir, 'unit_test_sample_coordinates.txt');
      data_dir_path = fullfile(work_dir, 'data');
      gibbs_data_file_name = 'sample_gibbs_data.mat';
      em_data_file_name = 'sample_em_ata.mat';
      system(['mkdir -p ' data_dir_path]);
      CBIG_AuthorTopicEM_GenerateDataFromText(text_data_path, ...
        data_dir_path, gibbs_data_file_name, em_data_file_name, ...
        brain_mask1mm, brain_mask2mm);
    
      data = load(fullfile(data_dir_path, em_data_file_name));
      assert(size(data.paradigm_by_exp, 1) == 7, 'Wrong number of paradigms');
      
      % inference
      all_Ks = 2:2; % set all_Ks to higher number of components to discover the full range of solutions
      alpha = 100;
      eta = 0.01;
      seeds = 1:3; % set seeds = 1:1000 to obtain stable estimates
    
      % gibbs sampler's parameters
      gibbs_burn_in = 2; % set tto 100 or more in actual experiments
      gibbs_inverval = gibbs_burn_in;
      num_gibbs_samples = 1; % set to 100 or more in actual experiments
      gibbs_input_path = fullfile(data_dir_path, gibbs_data_file_name);
    
      disp('Start Gibbs sampling');
      for K = all_Ks
        disp(['  K: ' num2str(K)]);
        CBIG_AuthorTopicEM_RunGibbsSampler(seeds, K, gibbs_burn_in, gibbs_inverval, num_gibbs_samples, ...
            gibbs_input_path, work_dir, brain_mask2mm, alpha, eta);
      end
      disp('Finished Gibbs sampling');
    
      em_init_type = 'GIBBS';
      em_input_path = fullfile(data_dir_path, em_data_file_name);
      disp('Start EM inference');
      for K = all_Ks
        for seed = seeds
          gibbs_init_dir = fullfile(work_dir, 'GIBBS_outputs', ['K' num2str(K)], ...
            ['alpha' num2str(alpha) '_eta' num2str(eta)]);
          gibbs_init_file = fullfile(gibbs_init_dir, ['seed' num2str(seed)], ...
            ['seed' num2str(seed) '.burn' num2str(gibbs_burn_in) '.int' ...
            num2str(gibbs_inverval) '.burnin.mat']);
          CBIG_AuthorTopicEM_RunEM(seed, K, em_input_path, brain_mask2mm, alpha, eta, em_init_type, gibbs_init_file);
        end
      end
      disp('Finished EM inference');
    
      % find the best solutions from all reinitializations
      output_dir = fullfile(work_dir, ['EM_outputs_' em_init_type '_init']);
      best_solution_dir = fullfile(output_dir, 'best_solution', ['alpha' num2str(alpha) '_eta' num2str(eta)]);
      for K = all_Ks
        CBIG_AuthorTopicEM_ComputeBestSolution(output_dir, K, seeds, alpha, eta);
    
        best_solution = load(fullfile(best_solution_dir, ['best_solution_K' num2str(K, '%03d') '.mat']));
        assert(size(best_solution.params.beta, 1) == K, 'Wrong number of components in beta');
        assert(size(best_solution.params.beta, 2) == sum(brain_mask2mm.vol(:) ~= 0), ...
          'Mismatched number of voxels in beta');
        assert(size(best_solution.params.theta, 2) == K, 'Wrong number of components in theta');
        assert(size(best_solution.params.theta, 1) == size(data.paradigm_by_exp, 1), ...
          'Mismatched number of tasks in theta');
      end
    
      % find the most typical solutions from all reinitializations
      avg_solution_dir = fullfile(output_dir, 'avg_solution', ['alpha' num2str(alpha) '_eta' num2str(eta)]);
      for K = all_Ks
        CBIG_AuthorTopicEM_ComputeMostTypicalSolution(output_dir, K, seeds, alpha, eta);
      
        avg_solution = load(fullfile(avg_solution_dir, ['avg_solution_K' num2str(K, '%03d') '.mat']));
        assert(size(avg_solution.params.beta, 1) == K, 'Wrong number of components in beta');
        assert(size(avg_solution.params.beta, 2) == sum(brain_mask2mm.vol(:) ~= 0), ...
          'Mismatched number of voxels in beta');
        assert(size(avg_solution.params.theta, 2) == K, 'Wrong number of components in theta');
        assert(size(avg_solution.params.theta, 1) == size(data.paradigm_by_exp, 1), ...
          'Mismatched number of tasks in theta')
      end
    
      % visualize components, e.g. 2-component best solution
      figures_dir = fullfile(work_dir, 'figures');
      input_path = fullfile(best_solution_dir, 'best_solution_K002.mat');
      CBIG_AuthorTopicEM_VisualizeComponentsOnBrainSurface(input_path, figures_dir);
      assert(exist(fullfile(figures_dir, 'clear_brain_min1e-5_max5e-5', ...
        'C1.grid.png'), 'file') == 2, 'Missing visualization for component C1');
      assert(exist(fullfile(figures_dir, 'clear_brain_min1e-5_max5e-5', ...
        'C2.grid.png'), 'file') == 2, 'Missing visualization for component C2');
    
      % clean up
      rmpath(fullfile(CBIG_CODE_DIR, 'external_packages', 'matlab', ...
        'non_default_packages', 'topictoolbox')); % Gibbs sampler
      rmpath(fullfile(utilities_dir, 'preprocessing'));
      rmpath(fullfile(utilities_dir, 'inference'));
      rmpath(fullfile(utilities_dir, 'visualization'));
      
      system(['rm -r ', data_dir_path]);
      system(['rm -r ', fullfile(work_dir, 'GIBBS_outputs')]);
      system(['rm -r ', output_dir]);
      system(['rm -r ', figures_dir]);
    end
  end
end
