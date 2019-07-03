function params = CBIG_AuthorTopic_EM_InitParams(exp_acts, paradigm_by_exp, params)
% params = CBIG_AuthorTopic_EM_InitParams(exp_acts, paradigm_by_exp, params)
%
% Initialize parameters of the author-topic model
%
% Input:
%   - exp_acts        = D x 1 cell where D is the number of experiments
%                       exp_acts{d} is a Nd x V sparse matrix, where Nd is the number of activation foci
%                       in the experiment, and V is the number of brain voxels. exp_acts{d}(n, j) = 1
%                       if nth activation focus of experiment d is j-th brain voxel.
%   - paradigm_by_exp = A x D logical sparse matrix, where A is the number of paradigms,
%                       D is the number of experiments.
%                       paradigm_by_exp(a, d) = 1 if paradigm "a" is in experiment d
%   - params          = struct of models' parameters declared in CBIG_AuthorTopic_EM_CreateEmptyParams
%
% Output:
%   - params          = initialized models' parameters
%
% Example:
%   params = CBIG_AuthorTopic_EM_InitParams(exp_acts, paradigm_by_exp, params)
%   Initialize parameters of the author-topic model with the activation foci of the experiments defined in
%   exp_acts, task paradigms of the experiments defined in paradigm_by_exp, and the initial parameters
%   stored in params.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  params.T = size(paradigm_by_exp, 1);
  params.V = size(exp_acts{1}, 2);
  params.E = length(exp_acts);
  
  if(sum(params.brain_mask.vol(:) == 1) ~= params.V)
     error('Brain mask is different size from dictionary'); 
  end
  
  if(strcmp(params.init_type, 'RAND'))
     disp('Initializing randomly');
     
     params.theta     = rand([params.T params.K]);
     params.theta     = params.theta;
     params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
     
     params.beta = rand([params.K params.V]);
     params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
  elseif(strcmp(params.init_type, 'GIBBS'))
     disp('Initializing using Gibbs input');
     load(params.init_file);
     
     % initialize theta
     if(size(theta, 1) ~= params.T || size(theta, 2) ~= params.K)
         error('Initialization theta not the same size');
     end
     params.theta = full(theta);
     
     % initialize beta
     for i = 1:params.K
         for j = 1:params.init_smooth
             process_maps.vol(:, :, :, i) = CBIG_Smooth3DVolumeWithMasks(squeeze(process_maps.vol(:, :, :, i)), ...
                                                                params.brain_mask.vol, 'SAME', 'box', 3);
         end
     end
     
     process_maps = reshape(process_maps.vol, ...
       [size(process_maps.vol, 1)*size(process_maps.vol, 2)*size(process_maps.vol, 3) ...
        size(process_maps.vol, 4)]);
     process_maps = process_maps';
     process_maps = process_maps(:, params.brain_mask.vol(:) == 1);
     
     if(size(process_maps, 1) ~= params.K || size(process_maps, 2) ~= params.V)
         error('Initialization beta not the same size');
     end
     params.beta  = full(process_maps);
     params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
  else
     error('No other initialization');
  end
  
  params.log_theta = log(params.theta); 
  params.log_beta = log(params.beta);
