function params = CBIG_AuthorTopicEM_Inference(act, paradigm_by_exp, params)
% params = CBIG_AuthorTopicEM_Inference(act, paradigm_by_exp, params)
%
% Runs Expectation-Maximization (EM) algorithm to estimate parameters of the author-topic model
%
% Input:
%   - act          = E x 1 cell where E is the number of experiments
%                     act{e} is a Ne x V sparse matrix, where Ne is the number of unique activation
%                     foci in experiment e, and V is the number of voxels. act{e}(n, j) = 1
%   - paradigm_by_exp = T x E logical sparse matrix, where T is the number of paradigms, and
%                     E is the number of experiments.
%                     paradigm_by_exp(t, e) = 1 if paradigm t is utilized in experiment e
%   - params       = struct of parameters used by the model
%                  params is created by CBIG_AuthorTopicEM_CreateEmptyParams
% Output:
%   - params       = updated struct with the model parameters estimated by EM algorithm
%
% Example:
%   params = CBIG_AuthorTopicEM_Inference(act, paradigm_by_exp, params)
%   Perform inference with the EM algorithm to estimate parameters
%   of the author-topic model and update the params struct.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  rand('twister', sum(10000*params.seed));
  
  params = CBIG_AuthorTopicEM_InitParams(act, paradigm_by_exp, params);
  paradigm_by_exp = full(paradigm_by_exp);
  
  params.log_likelihood = -inf;
  tic;
  for iter = 1:params.em_max_iter
      
      params.doc_log_likelihood = zeros(params.E, 1);
      params.new_theta = zeros([params.T params.K]);
      params.new_beta  = zeros([params.K params.V]);
      
      % e-step
      for e = 1:params.E
          if(size(act, 2) == 3)
              [params, params.doc_log_likelihood(e)] = CBIG_AuthorTopicEM_ExpEstep( ...
                act(e, :), paradigm_by_exp(:, e), params);
          else
              error('dimensions of act wrong');
          end
      end
      
      % m-step
      params.theta     = params.new_theta + params.virtual_theta;
      params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
      params.log_theta = log(params.theta);
      
      params.beta     = params.new_beta + params.virtual_beta;
      params.beta     = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
      params.log_beta = log(params.beta);
  
      % compute convergence conditions
      new_log_likelihood = sum(params.doc_log_likelihood);
      if(isinf(params.log_likelihood))
          log_likelihood_improvement = inf;
      else
          log_likelihood_improvement = abs(new_log_likelihood - params.log_likelihood)/abs(params.log_likelihood);
      end
      params.log_likelihood = new_log_likelihood;
      disp(['Iter ' num2str(iter, '%03d') ': log likelihood = ' ...
        num2str(params.log_likelihood) ' (' num2str(log_likelihood_improvement) '), time elapsed = ' num2str(toc)]);
      tic;
      
      if(log_likelihood_improvement < params.em_convergence && ...
         iter >= params.em_min_iter && mod(iter, params.mod_smooth) ~= 1)
          break;
      end
      
      if(mod(iter, params.mod_smooth) == 0 && iter ~= params.em_max_iter)
          disp('Smoothing ...');
          
          % reshape
          process_maps = zeros([ ...
            size(params.brain_mask.vol, 1)*size(params.brain_mask.vol, 2)*size(params.brain_mask.vol, 3) params.K]);
          for k = 1:params.K
              process_maps(params.brain_mask.vol(:) == 1, k) = params.beta(k, :);
          end
          process_maps = reshape(process_maps, [size(params.brain_mask.vol) params.K]);
          
          % smooth
          for k = 1:params.K
              for j = 1:params.init_smooth
                  process_maps(:, :, :, k) = CBIG_Smooth3DVolumeWithMasks(squeeze(process_maps(:, :, :, k)), ...
                                                                     params.brain_mask.vol, 'SAME', 'box', 3);
              end
          end
          
          process_maps = reshape(process_maps, ...
            [size(process_maps, 1)*size(process_maps, 2)*size(process_maps, 3) size(process_maps, 4)]);
          process_maps = process_maps';
          process_maps = process_maps(:, params.brain_mask.vol(:) == 1);
          
          params.beta  = full(process_maps);
          params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
          params.log_beta = log(params.beta);
      end
  end
  
  % final e-step to compute final log likelihood
  params.doc_log_likelihood = zeros(params.E, 1);
  for e = 1:params.E
      if(size(act, 2) == 3)
          q = CBIG_AuthorTopicEM_ExpInference(act(e, :), paradigm_by_exp(:, e), params);
          params.doc_log_likelihood(e) = CBIG_AuthorTopicEM_ExpLogLikelihood( ...
            act(e, :), paradigm_by_exp(:, e), params, q);
      else
          error('Dimensions of act variable is wrong');
      end
  end
  params.log_likelihood = sum(params.doc_log_likelihood);
  disp(['Final log likelihood: ' num2str(params.log_likelihood)]);
  
  params = rmfield(params, 'log_theta');
  params = rmfield(params, 'log_beta');
  params = rmfield(params, 'new_theta');
  params = rmfield(params, 'new_beta');
