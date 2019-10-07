function m = CBIG_nanstd(x, norm, dim)

% m = CBIG_nanstd(x, norm, dim)
% 
% Compute the std values of a matrix/vector excluding NaN entries.
% It has the same functionality as nanstd.m, but avoid to use statistics
% toolbox.
% 
% - Inputs
%   x   :    the data matrix/vector you want to compute the mean.
%   norm:    the norm is a binary scalar {0,1}. If norm = 0, this means
%            that the empirical std will be normalised by N-1 where N is
%            the number of observations. If norm = 1, this means that the
%            emprical std will be normalised by N.
%   dim :    takes the std along dimension dim of x. If dim is not
%            passed in, the std will be taken along the first 
%            non-singleton dimension. 
% 
% - Outputs
%   m   :    the output matrix/vector/scalar of the std excluding NaNs.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(nargin < 3)
    dim = 1;
    if(size(x, 1) == 1)
        dim = find(size(x) > 1, 1);
    end
end

% find nan entries
nanind = isnan(x);

% set nan entries to be zeros and sum
xmean = CBIG_nanmean(x, dim);
x_norm = bsxfun(@minus,x, xmean);
x_norm(nanind) = 0;
x_sq = x_norm.^2;
xsum = sum(x_sq, dim);

% count nan-nan entries
if norm == 0
    count = size(x, dim) - sum(nanind, dim) - 1;
else 
    count = size(x, dim) - sum(nanind, dim);
end

% nanmean
m = xsum ./ count;
m = m.^(0.5);


end