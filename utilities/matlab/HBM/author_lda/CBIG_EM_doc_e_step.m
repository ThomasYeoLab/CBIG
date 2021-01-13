function [params, doc_log_likelihood, q] = CBIG_EM_doc_e_step(w, paradigm, params) 

% [params, doc_log_likelihood, q] = CBIG_EM_doc_e_step(w, paradigm, params)
%
% Perform E-step of the Author-Topic model on each experiment (document). Note that
% compared to CBIG_EM_doc_e_step_wc.m, this function takes in an argument w of a different structure
% FORMAT [params, doc_log_likelihood, q] = CBIG_EM_doc_e_step(w, paradigm, params)
%
% w        = Nd x V sparse matrix, where Nd is the number of unique activation foci
%            (unique words) in the document and V is the number of voxels (vocabulary
%            words). w(n, v) = 1 if the n-th activatin foci (word) in the experiment
%            (document) is the v-th voxel (vocabulary word)
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

if size(paradigm,2) ~= 1
    error('Input argument ''paradigm'' should be a column vector');
end

q = CBIG_EM_doc_inference(w, paradigm, params);

% Compute document log likelihood
doc_log_likelihood = CBIG_EM_doc_log_likelihood(w, paradigm, params, q);

% update theta for M-step: theta is A x T
params.new_theta(paradigm, :) = params.new_theta(paradigm, :) + squeeze(sum(q, 3)); % Ad x T

% update beta for M-step: beta is T x V
beta_update = squeeze(sum(q, 1)); % T x Nd
params.new_beta = params.new_beta + beta_update * w;

