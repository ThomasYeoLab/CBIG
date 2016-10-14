function div = CBIG_KLdiv(ref_prob, prob)

% Calculate KL divergence.
%
% 	div = CBIG_KLdiv(ref_prob, prob)
% 	Input: 
% 		ref_prob: reference probability, Nx1 vector
% 		prob	: probability, Nx1 vector
% 	Output:
%		div		: KL divergence
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


div = sum(ref_prob .* (log2(ref_prob) - log2(prob)));
