function varargout = palm_boxcox(varargin)
% Find the lambda parameter and compute the Box-Cox transformation.
%
% Usage:
% [Y,L] = palm_boxcox(X)
% Y     = palm_boxcox(X,L)
% 
% X : Input data. All values need to be positive.
% Y : Transformed data.
% L : Lambda parameter.
%
% Reference:
% Box GEP and Cox DR. An Analysis of Transformations.
% J R Stat Soc Series B. 1964;26(2):211-252.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

X = varargin{1};
if nargin == 1,
    Lfun = @(L)loglike(X,L);
    L    = fminsearch(Lfun,0,optimset('Display','off','MaxFunEvals',1000));
    varargout{1} = bct(X,L);
else
    L = varargin{2};
    varargout{1} = bct(X,L);
end
if nargout == 2,
    varargout{2} = L;
end

% ==============================================================
function LL = loglike(X,L)
% Likelihood function.
Y  = bct(X,L);
LL = numel(X)*log(std(Y,1,1)^2)/2 - (L-1)*sum(log(X));

% ==============================================================
function Y = bct(X,L)
% Actual Box-Cox transformation.
if L == 0,
    Y = log(X);
else
    Y = (X.^L - 1)/L;
end
