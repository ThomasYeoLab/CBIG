function params = CBIG_AuthorTopicEM_CreateEmptyParams(seed, K, brain_mask2mm, init_type, init_file)
% params = CBIG_AuthorTopicEM_CreateEmptyParams(seed, K, brain_mask2mm, init_type, init_file)
%
% Generate a struct of parameters for the author-topic model
%
% Input:
%   - seed          = seed used by Matlab's random number generator
%   - K             = number of components to be estimated
%   - brain_mask2mm = struct containing a brain mask at 2mm resolution (read by MRIread)
%   - init_type     = a string used to determine how the auxiliary parameter
%                     of the model is generated. 'RAND' if the auxiliary
%                     parameter is randomly generated, 'GIBBS' if the auxiliary
%                     parameter is the output of Gibbs sampling
%   - init_file     = if init_type = 'GIBBS', the output of GIBBS sampling from
%                     init_file is used
% Output:
%   - params        = struct of parameters used by the Author-Topic model
%
% Example:
%   init_file = '~/Work/GIBBS_outputs/K2/alpha100_eta0.01/seed1/seed1.burn2.int2.sample1.mat';
%   work_dir = '~/Work';
%   brain_mask2mm = MRIread(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
%     'meta-analysis', 'Yeo2015_AuthorTopicEM', 'utilities', ...
%     'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz'));
%   CBIG_AuthorTopicEM_CreateEmptyParams(42, 2, brain_mask2mm, 'GIBBS', init_file)
%  
%   Set up parameters for the author-topic model with 2 components, initialization seed 42
%   a 2mm brain mask read by MRIread. The Expectation-Maximization (EM) algorithm is initialized
%   with output from Gibbs sampler.
%  
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  params.brain_mask     = brain_mask2mm;
  params.em_max_iter    = 1000;
  params.em_min_iter    = 10;
  params.em_convergence = 1e-5;
  params.num_smooth     = 1;
  params.mod_smooth     = 1000;
  params.max_theta      = 1.1;
  
  if(nargin < 4)
      params.init_type = 'RAND';
  else
      params.init_type = init_type;
      params.init_file = init_file;
      params.init_smooth = 5;
  end
  
  params.seed = seed;
  params.K = K;
