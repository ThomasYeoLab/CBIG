function q = CBIG_EM_doc_inference(w, paradigm, params) 

% q = CBIG_EM_doc_inference(w, paradigm, params) 
%
% Update the auxillary parameter q in the E-step on each experiment (document).
% Note that compared to CBIG_EM_doc_inference_wc.m, this function takes in an argument
% w of a different structure
% FORMAT q = CBIG_EM_doc_inference(w, paradigm, params)
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
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(paradigm,2) ~= 1
    error('Input argument ''paradigm'' should be a column vector');
end

num_paradigms = sum(paradigm);

theta = params.theta(paradigm, :); % num_paradigms x T
beta_wd = params.beta * w'; % T x Nd
q = zeros([num_paradigms params.T size(w, 1)]);
for a = 1:num_paradigms
   q(a, :, :) = bsxfun(@times, beta_wd, theta(a, :)'); 
end
%q = bsxfun(@times, repmat(reshape(beta_wd, [1 params.T size(w, 1)]), [num_paradigms 1 1]), theta); 

normalizer = 1./squeeze(sum(squeeze(sum(q, 2)), 1));
%q = bsxfun(@times, q, reshape(normalizer, [1 1 length(normalizer)]));
for a = 1:num_paradigms
   q(a, :, :) = bsxfun(@times, squeeze(q(a, :, :)), normalizer); 
end
