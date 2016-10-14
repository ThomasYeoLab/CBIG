function P = CBIG_chi2Tests_PE(X)

% Computes p value of pearson chi square tests (fully vectorized format).
% 
%   P = CBIG_chi2Tests_PE(X)
%   Input: 
%       X  : M_1 x M_2 x N matrix
%   Output: 
%       P  : 1 x N vector, p value
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


rs = sum(X, 1);
cs = sum(X, 2);

Ns = sum(rs, 2);

rep_rs = repmat(rs, [size(X, 1) 1 1]);
rep_cs = repmat(cs, [1 size(X, 2) 1]);

Eo = bsxfun(@times, rep_rs .* rep_cs, 1./Ns);

t = squeeze(sum(sum((Eo-X).^2 ./ Eo, 2), 1));

df=prod([size(X,1) size(X,2)]-[1,1]); % degree of freedom
P = 1-chi2cdf(t,df);

