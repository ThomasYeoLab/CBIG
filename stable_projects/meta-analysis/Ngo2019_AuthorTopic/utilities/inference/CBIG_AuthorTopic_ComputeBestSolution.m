function bestSeed = CBIG_AuthorTopic_ComputeBestSolution(rootDir, K, ...
    seeds, alpha, eta)
% bestSeed = CBIG_AuthorTopic_ComputeBestSolution(rootDir, K, ...
%   seeds, alpha, eta)
%
% Return the index of the initialization (seed) of the Collapsed
% Variational Bayes (CVB) algorithm run that returned the highest
% variational bound of the data log likelihood of the author-topic model
% parameters' estimates.
%
% Input:
%  - rootDir    : base directory containing all the estimates of the
%                 author-topic model parameters by the CVB algorithm across
%                 all initializations (seeds) at K components
%                 and with hyperparameters alpha and eta of the Dirichlet
%                 priors. The paths to the model estimates have the
%                 following format:
%  <rootDir>/K<K>/alpha<alpha>_eta<eta>/params_K<K>_SEED<seed number>.mat
%  - K: number of cognitive components/ co-activation patterns.
%  - seeds      : array of numbers denoting the indices of the
%                 initializations (seeds).
%  - alpha & eta: hyperparameters of the author topic model's Dirichlet
%                 priors.
% Output:
%  - bestSeed: index of the initialization (seed) that returns the highest
%              variational bound of the data log likelihood.
%              A copy of the best parameters estimate is also saved at
%    <rootDir>/outputs/bestSolution/BestSolution_K<K>.mat
%
% Example:
%   bestSeed = CBIG_AuthorTopic_ComputeBestSolution('/Work/outputs', 2, ...
%              [1:100], 100, 0.01)
%   Return the initialization (seed) that produces the highest variational
%   bound of the data log likelihood of the author-topic model's parameter
%   estimates. The indices of the candidate seeds are from 1 to 100,
%   inclusive. The number of cognitive components/ co-activation patterns
%   is 2, and the model's hyperparameters are 100 and 0.01.
%   A copy of the best estimate is saved at
%   '/Work/ouputs/bestSolution/BestSolution_K002.mat'
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % load mask
  alphaEtaPrefix = ['alpha' num2str(alpha) '_eta' num2str(eta)];
  inputDir = fullfile(rootDir, ['K' num2str(K)], alphaEtaPrefix);

  varBound = zeros(length(seeds), 1);

  % load samples
  disp('loading samples ...');
  count = 1;
  for seed = seeds
      disp(['seed: ' num2str(seed)]);
      load(fullfile(inputDir, ['params_K' num2str(K) '_SEED' num2str(seed) '.mat']));
      mask = MRIread(params.maskPath);
      maskIndices = find(mask.vol > 0);
      if numel(maskIndices) ~= size(params.beta, 2)
        error('Brain mask is different size from beta');
      end;

      if(max(params.theta(:)) >= 1.1)
          error('Invalid theta estimates');
      else
          varBound(count) = params.varBound(numel(params.varBound));
      end
      count = count + 1;
  end

  [~, bestSolution] = max(varBound);

  outputDir = fullfile(rootDir, 'bestSolution', ['alpha' num2str(params.alpha) '_eta' num2str(params.eta)]);
  mkdir(outputDir);

  % save best solution
  seed = seeds(bestSolution);
  load(fullfile(inputDir, ['params_K' num2str(K) '_SEED' num2str(seed) '.mat']));
  save(fullfile(outputDir, ['BestSolution_K' num2str(K, '%03d') '.mat']), 'params');
  bestSeed = seed;
