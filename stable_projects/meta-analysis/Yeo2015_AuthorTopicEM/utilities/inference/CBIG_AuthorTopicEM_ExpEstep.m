function [params, exp_log_likelihood, q] = CBIG_AuthorTopicEM_ExpEstep(w, paradigm, params) 
% [params, exp_log_likelihood, q] = CBIG_AuthorTopicEM_ExpEstep(w, paradigm, params) 
%
% Perform E-step of the Author-Topic model on a given experiment.
%
% Input:
%   - w           = 1 x 3 cell array of experiment's activation foci.
%     w{1}        = Nd x V sparse matrix where Nd is the number of unique activation foci
%                   in the given experiment, and V is the number of brain voxels.
%                   w{1} is the same as w{2} but has counts.
%     w{2}        = Nd x V sparse matrix where w{2}(n, v) = 1 if the n-th activation focus in
%                   the experiment is the v-th voxel.
%     w{3}        = 1 x Nd vector, where w{n} is the number of times the n-th unique activation
%                 in the experiment occurs.
%  Note that w{1} = bsxfun(@times, w{2}, w{3}');
%   - paradigm    = T x 1 logical sparse matrix, where T is the number of task paradigms.
%                 paradigm(t) = 1 if task "t" is recruited in the given experiment.
%   - params      = struct of model's parameters
%
% Output:
%   - params      = updated model's parameters.
%   - exp_log_likelihood = log likelihood of a single experiment.
%   - q           = Te x K x Nd matrix defining the auxiliary parameter of the EM algorithm,
%                where Te is the number of task paradigms of experiment e,
%                K is the number of components of the author-topic model,
%                and Ne is the number of unique activation foci reported in experiment e.
%
% Example:
%   [params, exp_log_likelihood, q] = CBIG_AuthorTopic_EM_ExpEstep(w, paradigm, params) 
%   Perform E-step update on an individual experiment with the activation pattern defined by
%   w, the experiment's task defined by paradigm and the model's parameters params. Returned
%   variable is the updated params variable, log likelihood exp_log_likelihood
%   and the updated auxiliary parameter q specific to the given experiment.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  q = CBIG_AuthorTopicEM_ExpInference(w, paradigm, params);
  
  % Compute document log likelihood
  exp_log_likelihood = CBIG_AuthorTopicEM_ExpLogLikelihood(w, paradigm, params, q);
  
  % update theta for M-step: theta is T x K
  params.new_theta(paradigm, :) = params.new_theta(paradigm, :) + ...
    squeeze(sum(bsxfun(@times, q, reshape(w{3}, [1 1 length(w{3})])), 3));
  
  % update beta for M-step: beta is K x V
  beta_update = squeeze(sum(q, 1)); % K x Ne
  
  params.new_beta = params.new_beta + beta_update * w{1};
