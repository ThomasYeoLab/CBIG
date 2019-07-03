function q = CBIG_AuthorTopic_EM_Exp_Inference(w, paradigm, params) 
% q = CBIG_AuthorTopic_EM_Exp_Inference(w, paradigm, params) 
%
% Update auxiliary parameter q in the E-step on each experiment
%
% Input:
%   - w           = 1 x 3 cell array of experiment's activation foci.
%     w{1}        = Ne x V sparse matrix where Ne is the number of unique activation foci
%                   in the given experiment, and V is the number of brain voxels.
%                   w{1} is the same as w{2} but has counts.
%     w{2}        = Ne x V sparse matrix where w{2}(n, v) = 1 if the n-th activation focus in
%                   the experiment is the v-th voxel.
%     w{3}        = 1 x Ne vector, where w{n} is the number of times the n-th unique activation
%                 in the experiment occurs.
%  Note that w{1} = bsxfun(@times, w{2}, w{3}');
%   - paradigm    = T x 1 logical sparse matrix, where T is the number of task paradigms.
%                 paradigm(t) = 1 if task "t" is recruited in the given experiment.
%   - params      = struct of model's parameters
%
% Output:
%   - q        = Te x K x Ne matrix
%                where K is the number of components of the author-topic model,
%                and Te and Ne are the number of task paradigms and number of unique activation
%                foci reported in the current experiment e respectively.
%
% Example:
%   q = CBIG_AuthorTopic_EM_Exp_Inference(w, paradigm, params)
%   Return updated auxiliary parameter q based on a given experiment with the activation pattern
%   defined by w, the experiment's task defined by paradigm and the model's parameters in params.
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
  num_paradigms = sum(paradigm);
  
  theta = params.theta(paradigm, :); % num_paradigms x K
  beta_wd = params.beta * w{2}'; % (K x Nd)*
  q = zeros([num_paradigms params.K length(w{3})]);
  for t = 1:num_paradigms
     q(t, :, :) = bsxfun(@times, beta_wd, theta(t, :)'); 
  end
  
  normalizer = 1./squeeze(sum(squeeze(sum(q, 2)), 1));
  for t = 1:num_paradigms
     q(t, :, :) = bsxfun(@times, squeeze(q(t, :, :)), normalizer); 
  end
