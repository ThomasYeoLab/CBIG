function varargout = palm_yeojohnson(varargin)
% Find the lambda parameter and compute the Yeo-Johnson transformation.
%
% Usage:
% [Y,L] = palm_yeojohnson(X)
% Y     = palm_yeojohnson(X,L)
% 
% X : Input data.
% Y : Transformed data.
% L : Lambda parameter.
%
% Reference:
% Yeo IK and Johnson RA. A new family of power transformations to
% improve normality or symmetry. Biometrika. 2000;87(4):954-9.
% _____________________________________
% Anderson M. Winkler
% Jul/2017
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
    L    = fminsearch(Lfun,1,optimset('Display','off','MaxFunEvals',1000));
    varargout{1} = yjt(X,L);
else
    L = varargin{2};
    varargout{1} = yjt(X,L);
end
if nargout == 2,
    varargout{2} = L;
end

% ==============================================================
function LL = loglike(X,L)
% Log-likelihood function (eq 3.1 of the paper), with the
% sign reversed to be minimised.
Y   = yjt(X,L);
N   = numel(X);
mu  = mean(Y);
sg2 = var(Y,1);
LL  = N*(log(2*pi*sg2))/2   ...
    + sum((Y-mu).^2)/2/sg2  ...
    - (L-1)*sum(sign(X.*log(abs(X)+1)));

% ==============================================================
function Y = yjt(X,L)
% Yeo-Johnson transformation (eq 2.1 of the paper)
Y   = nan(size(X));
pos = X >= 0;
if L == 0,
    Y( pos) = log(X(pos)+1);                  % X >= 0, L == 0
    Y(~pos) = -((1-X(~pos)).^(2-L)-1)/(2-L);  % X <  0, L ~= 2
elseif L == 2,
    Y( pos) = ((X(pos)+1).^L-1)./L;           % X >= 0, L ~= 0
    Y(~pos) = -log(1-X(~pos));                % X <  0, L == 2
else
    Y( pos) = ((X(pos)+1).^L-1)./L;           % X >= 0, L ~= 0
    Y(~pos) = -((1-X(~pos)).^(2-L)-1)/(2-L);  % X <  0, L ~= 2
end
