function m = CBIG_nanmean(x, dim)

% m = CBIG_nanmean(x, dim)
% 
% Compute the mean values of a matrix/vector excluding NaN entries.
% It has the same functionality as nanmean.m, but avoid to use statistics
% toolbox.
% 
% - Inputs
%   x   :    the data matrix/vector you want to compute the mean.
%   dim :    takes the mean along dimension dim of x. If dim is not
%            passed in, the mean will be taken along the first 
%            non-singleton dimension. 
% 
% - Outputs
%   m   :    the output matrix/vector/scalar of the mean excluding NaNs.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(nargin == 1)
    dim = 1;
    if(size(x, 1) == 1)
        dim = find(size(x) > 1, 1);
    end
end

% find nan entries
nanind = isnan(x);

% set nan entries to be zeros and sum
x(nanind) = 0;
xsum = sum(x, dim);

% count nan-nan entries
count = size(x, dim) - sum(nanind, dim);

% nanmean
m = xsum ./ count;


end