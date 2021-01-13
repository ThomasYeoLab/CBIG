function q = CBIG_EM_doc_inference_wc(w, paradigm, params) 

% q = CBIG_EM_doc_inference_wc(w, paradigm, params) 
%
% Update auxillary parameter q in the E-step on each experiment (document). Note that
% compared to CBIG_EM_doc_inference.m, this function takes in an argument w of a different structure
% FORMAT q = CBIG_EM_doc_inference_wc(w, paradigm, params)
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
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(w{3},1) ~= 1
    error('Input argument ''w{3}'' should be a row vector');
end
if size(paradigm,2) ~= 1
    error('Input argument ''paradigm'' should be a column vector');
end

num_paradigms = sum(paradigm);

theta = params.theta(paradigm, :); % num_paradigms x T
beta_wd = params.beta * w{2}'; % T x Nd*
q = zeros([num_paradigms params.T length(w{3})]);
for a = 1:num_paradigms
   q(a, :, :) = bsxfun(@times, beta_wd, theta(a, :)'); 
end
%q = bsxfun(@times, repmat(reshape(beta_wd, [1 params.T size(w, 1)]), [num_paradigms 1 1]), theta); 

normalizer = 1./squeeze(sum(squeeze(sum(q, 2)), 1));
%q = bsxfun(@times, q, reshape(normalizer, [1 1 length(normalizer)]));
for a = 1:num_paradigms
   q(a, :, :) = bsxfun(@times, squeeze(q(a, :, :)), normalizer); 
end
