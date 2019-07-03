function best_seed = CBIG_AuthorTopicEM_ComputeBestSolution(root_dir, K, ...
    seeds, alpha, eta)
% best_seed = CBIG_AuthorTopicEM_ComputeBestSolution(root_dir, K, ...
%   seeds, alpha, eta)
%
% Return the index of the initialization (seed) of the Expectation-
% Maximizaiton (EM) inference algorithm run that returned the highest
% data log likelihood given the author-topic model parameters' estimates.
%
% Input:
%  - root_dir    = base directory containing all the estimates of the
%                 author-topic model parameters by the EM algorithm across
%                 all initializations (seeds) at K components
%                 and with hyperparameters alpha and eta of the Dirichlet
%                 priors. The paths to the model estimates have the
%                 following format:
%  <root_dir>/K<K>/alpha<alpha>_eta<eta>/params_K<K>_SEED<seed number>.mat
%  - K           = number of cognitive components.
%  - seeds       = array of numbers denoting the indices of the
%                 initializations (seeds).
%  - alpha & eta = hyperparameters of the author topic model's Dirichlet
%                 priors.
% Output:
%  - best_seed =  index of the initialization (seed) that returns the highest
%              variational bound of the data log likelihood.
%              A copy of the best parameters estimate is also saved at
%    <root_dir>/outputs/best_solution/BestSolution_K<K>.mat
%
% Example:
%   best_seed = CBIG_AuthorTopicEM_ComputeBestSolution('/Work/outputs', 2, ...
%              [1:100], 100, 0.01)
%   Return the initialization (seed) that produces the highest variational
%   bound of the data log likelihood of the author-topic model's parameter
%   estimates. The indices of the candidate seeds are from 1 to 100,
%   inclusive. The number of cognitive components is 2, and the model's
%   hyperparameters are 100 and 0.01.
%   A copy of the best estimate is saved at
%   '/Work/outputs/best_solution/best_solution_K002.mat'
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  alpha_eta_prefix = ['alpha' num2str(alpha) '_eta' num2str(eta)];
  input_dir = fullfile(root_dir, ['K' num2str(K)], alpha_eta_prefix);

  likelihood = zeros(length(seeds), 1);

  % load samples
  disp('loading samples ...');
  count = 1;
  for seed = seeds
      disp(['seed: ' num2str(seed)]);
      load(fullfile(input_dir, ['seed' num2str(seed) '.mat']));
      mask_indices = find(params.brain_mask.vol > 0);

      assert(numel(mask_indices) == size(params.beta, 2), 'Brain mask has a different size from beta');
      assert(max(params.theta(:)) < 1.1, 'Invalid theta estimates');

      likelihood(count) = params.log_likelihood;
      count = count + 1;
  end

  [~, best_solution] = max(likelihood);

  output_dir = fullfile(root_dir, 'best_solution', alpha_eta_prefix);
  mkdir(output_dir);

  % save best solution
  seed = seeds(best_solution);
  load(fullfile(input_dir, ['seed' num2str(seed) '.mat']));
  save(fullfile(output_dir, ['best_solution_K' num2str(K, '%03d') '.mat']), 'params');
  best_seed = seed;
