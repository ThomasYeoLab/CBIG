function [params, doc_log_likelihood, q] = CBIG_EM_doc_e_step_wc(w, paradigm, params) 

% [params, doc_log_likelihood, q] = CBIG_EM_doc_e_step_wc(w, paradigm, params) 
%
% Perform E-step of the Author-Topic model on each experiment (document). Note that
% compared to CBIG_EM_doc_e_step.m, this function takes in an argument w of a different structure
% FORMAT [params, doc_log_likelihood, q] = CBIG_EM_doc_e_step_wc(w, paradigm, params)
%
% w           = 1 x 3 cell array
% w{1}        = Nd x V sparse matrix where Nd is the number of activation foci (unique
%               words) in the experiemnt (document), and V is the number of voxels
%               (vocabulary words). w{1} is the same as w{2} but has counts
% w{2}        = Nd x V sparse matrix where w{2}(n, v) = 1 if the n-th activation foci in i
%               the experiment is the v-th voxel
% w{3}        = 1 x Nd vector, where w{n} is the number of times the n-th unique activation
%               (word) in the experiment (document) appears
% Note that w{1} = bsxfun(@times, w{2}, w{3}');
% paradigm = A x 1 logical sparse matrix, where A is the number of paradigms (authors).
%            paradigm(a) = 1 if paradigm (author) "a" is in the current experiment (document)
% params   = struct of common parameters used by the model
%
% q        = Ad x T x Nd matrx, where Ad, and Nd are the number of paradigms (authors)
%            and unique activation foci (words) in the current experiment (document).
%            T is the number of components (topics)
% doc_log_likelihood = log likelihood of a single experiment (document)
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(w{3},1) ~= 1
    error('Input argument ''w{3}'' should be a row vector');
end
if size(paradigm,2) ~= 1
    error('Input argument ''paradigm'' should be a column vector');
end

q = CBIG_EM_doc_inference_wc(w, paradigm, params);

% Compute document log likelihood
doc_log_likelihood = CBIG_EM_doc_log_likelihood_wc(w, paradigm, params, q);

% update theta for M-step: theta is A x T
params.new_theta(paradigm, :) = params.new_theta(paradigm, :) ...
    + squeeze(sum(bsxfun(@times, q, reshape(w{3}, [1 1 length(w{3})])), 3)); % Ad x T

% update beta for M-step: beta is T x V
beta_update = squeeze(sum(q, 1)); % T x Nd

if(length(w{3}) == 1)
    params.new_beta = params.new_beta + beta_update' * w{1};
else
    params.new_beta = params.new_beta + beta_update * w{1};
end

