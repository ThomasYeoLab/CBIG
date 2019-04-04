function xsum = CBIG_nansum(x, dim)

% m = CBIG_nansum(x, dim)
% 
% Compute the sum of a matrix/vector excluding NaN entries.
% It has the same functionality as nansum.m, but avoid to use statistics
% toolbox.
% 
% - Inputs
%   x   :    the data matrix/vector you want to compute the sum.
%   dim :    takes the sum along dimension dim of x. If dim is not
%            passed in, the mean will be taken along the first 
%            non-singleton dimension. 
% 
% - Outputs
%   xsum   :    the output matrix/vector/scalar of the sum excluding NaNs.
% 
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

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




end