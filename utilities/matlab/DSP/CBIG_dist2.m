function out = CBIG_dist2(a,b)

% out = CBIG_dist2(a,b)
%
% For a: n x d, b: m x d, computes the distance matrix d_ij = ||a_i-b_j||^2
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



d = size(a,2);
out = repmat( (a.^2)*ones(d,1) ,1,size(b,1)) + repmat( (b.^2)*ones(d,1) ,1,size(a,1))' - 2 * a*b';