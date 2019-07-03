function exp_log_likelihood = CBIG_AuthorTopicEM_ExpLogLikelihood(w, paradigm, params, q) 
% exp_log_likelihood = CBIG_AuthorTopicEM_ExpLogLikelihood(w, paradigm, params, q) 
%
% Compute log likelihood of a given experiment given the auxiliary parameter of the EM algorithm
% and current model's parameters.
%
% Input:
%   - w           = 1 x 3 cell array of experiment's activation foci.
%     w{1}        = Ne x V sparse matrix where Ne is the number of unique activation foci
%                   in the given experiemnt, and V is the number of brain voxels.
%                   w{1} is the same as w{2} but has counts.
%     w{2}        = Nd x V sparse matrix where w{2}(n, v) = 1 if the n-th activation focus in
%                   the experiment is the v-th voxel.
%     w{3}        = 1 x Nd vector, where w{n} is the number of times the n-th unique activation
%                 in the experiment occurs.
%  Note that w{1} = bsxfun(@times, w{2}, w{3}');
%   - paradigm    = T x 1 logical sparse matrix, where T is the number of task paradigms.
%                 paradigm(t) = 1 if task "t" is recruited in the given experiment.
%   - params      = struct of model's parameters
%   - q           = updated auxiliary parameters of EM algorithm
%
% Output:
%   - exp_log_likelihood = log likelihood of the given experiment.
%
% Example:
%   exp_log_likelihood = CBIG_AuthorTopic_EM_Exp_LogLikelihood(w, paradigm, params, q)
%   Compute the log likelihood of an individual experiment with the activation pattern defined by
%   w, the experiment's task defined by paradigm, the model's parameters in params
%   and the auxiliary parameter q.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  log_likelihood_theta = sum(sum(sum(bsxfun(@times, bsxfun(@times, ...
    q, params.log_theta(paradigm, :)), reshape(w{3}, [1 1 length(w{3})])), 3), 2));
  
  beta_update = squeeze(sum(q, 1)); % (T x Ne)*
  log_likelihood_beta  = sum(sum(bsxfun(@times, ...
    beta_update .* (params.log_beta * w{2}'), reshape(w{3}, [1 length(w{3})])), 2));
  
  q_entropy = -sum(sum(sum(bsxfun(@times, q .* log(q), reshape(w{3}, [1 1 length(w{3})])), 3), 2));
  
  exp_log_likelihood = log_likelihood_theta + log_likelihood_beta + q_entropy; 
