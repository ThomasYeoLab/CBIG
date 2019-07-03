function avg_seed = CBIG_AuthorTopicEM_ComputeMostTypicalSolution(root_dir, K, seeds, alpha, eta)
% avg_seed = CBIG_AuthorTopicEM_ComputeMostTypicalSolution(root_dir, K, ...
%   seeds, alpha, eta)
%
% Return the index of the initialization (seed) of the Expectation-
% Maximization (EM) inference algorithm run that returned the author-topic
% model parameters' estimates that are closest (most typical) to the other
% initializations, measured in Pearson's correlation coefficient.
%
% Input:
%  - root_dir    : base directory containing all the estimates of the
%                 author-topic model parameters by the EM algorithm across
%                 all initializations (seeds) at K components
%                 and with hyperparameters alpha and eta of the Dirichlet
%                 priors. The paths to the model estimates have the
%                 following format:
%  <root_dir>/K<K>/alpha<alpha>_eta<eta>/params_K<K>_SEED<seed number>.mat
%  - K: number of cognitive components.
%  - seeds      : array of numbers denoting the indices of the
%                 initializations (seeds).
%  - alpha & eta: hyperparameters of the author topic model's Dirichlet
%                 priors.
% Output:
%  - avg_seed: index of the initialization (seed) whose estimates are the
%              closest to all other initializations.
%              A copy of the best parameters estimate is also saved at
%    <root_dir>/outputs/avg_solution/avg_solution_K<K>.mat
%
% Example:
%   best_seed = CBIG_AuthorTopicEM_ComputeMostTypicalSolution('/Work/outputs', 2, ...
%              [1:100], 100, 0.01)
%   Return the initialization (seed) that produces the author-topic model's parameter
%   estimates that are closest to those produced by other initializations. The indices
%   of the candidate seeds are from 1 to 100 inclusive. The number of cognitive
%   components is 2, and the model's hyperparameters are 100 and 0.01.
%   A copy of the most typical estimate is saved at
%   '/Work/ouputs/avg_solution/avg_solution_K002.mat'
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  alpha_eta_prefix = ['alpha' num2str(alpha) '_eta' num2str(eta)];
  input_dir = fullfile(root_dir, ['K' num2str(K)], alpha_eta_prefix);

  load(fullfile(input_dir, ['seed' num2str(seeds(1)) '.mat']));
  mask_indices = find(params.brain_mask.vol > 0);
  
  process_mat = zeros(length(mask_indices), K, length(seeds));
  
  % load samples
  disp('loading samples ...');
  count = 1;
  for seed = seeds
      disp(['seed: ' num2str(seed)]);
  
      load(fullfile(input_dir, ['seed' num2str(seed) '.mat']));
      mask_indices = find(params.brain_mask.vol > 0);
      if numel(mask_indices) ~= size(params.beta, 2)
        error('Brain mask is different size from beta');
      end;
  
      process_mat(:, :, count) = params.beta';
      count = count + 1;
  end
  
  disp('Computing correlation'); 
  process_mat_reshape = reshape(process_mat, [length(mask_indices) K*length(seeds)]);
  z = CBIG_self_corr(process_mat_reshape);
  
  disp('Computing assignments');
  cost = nan(length(seeds));
  for i = 1:length(seeds)
      for j = 1:length(seeds)
          if(i ~= j)
             row_start = (i - 1)*K + 1;
             row_end   = i*K;
             col_start = (j - 1)*K + 1;
             col_end   = j*K;
             x = z(row_start:row_end, col_start:col_end); 
             
             [tmp, cost(i, j)] = munkres(1 - x);  
          end
      end
  end
  avg_cost = nanmean(cost, 2);
  [Y, avg_seed] = min(avg_cost)
  
  % save avg solution
  output_dir = fullfile(root_dir, 'avg_solution', alpha_eta_prefix);
  mkdir(output_dir);
  
  load(fullfile(input_dir, ['seed' num2str(avg_seed) '.mat']));
  save(fullfile(output_dir, ['avg_solution_K' num2str(K, '%03d') '.mat']), 'params');
