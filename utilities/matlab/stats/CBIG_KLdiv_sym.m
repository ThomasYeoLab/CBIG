function div = CBIG_KLdiv_sym(prob1, prob2)

% Symmetric KL divergence, also called Jensen-Shannon divergence.
%
% 	div = CBIG_KLdiv_sym(prob1, prob2)
% 	Input: 
% 		prob1: 1st probability, Nx1 vector
% 		prob2: 2nd probability, Nx1 vector
% 	Output:
%		div		: symmetric KL divergence
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



div = (CBIG_KLdiv(prob1, prob2) + KLdiv(prob2, prob1))/2;
