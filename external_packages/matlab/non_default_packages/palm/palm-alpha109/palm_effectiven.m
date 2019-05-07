function n = palm_effectiven(nP,EE,ISE,islog)
% Given a number of shufflings, compute the effective
% sample size.
% 
% Usage:
% n = palm_effectiven(nP,EE,ISE,islog)
% 
% Inputs:
% - nP    : Number of shufflings (permutations, sign-flippings
%           or permutations with sign-flippings).
% - EE    : Boolean indicating whether the shufflings include
%           permutations (exchangeable errors).
%           Default is true.
% - ISE   : Boolean indicating whether the shufflings include
%           sign-flippings (independent and symmetric errors).
%           Default is false.
% - islog : Boolean indicating whether nP is given as log(nP).
%           Default is false.
% 
% Outputs:
% - n     : Effective sample size.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2014
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

% Simple argument parsing
narginchk(1,4);
if nargin == 1,
    EE    = true;
    ISE   = false;
    islog = false;
elseif nargin == 2,
    ISE   = false;
    islog = false;
elseif nargin == 3,
    islog = false;
end

% The Lambert's W function requires the 'specfun' package in Octave.
% For Matlab, it requires the Symbolic Math Toolbox.
if palm_isoctave,
    pkg load specfun
end

if EE && ~ISE,
    
    % Permutations only
    if islog,
        cte = nP - log(sqrt(2*pi));
        n   = cte./lambertw(cte/exp(1))-.5;
    else
        cte = log(nP/sqrt(2*pi));
        n   = cte./lambertw(cte/exp(1))-.5;
    end
    
elseif ~EE && ISE,
    
    % Sign-flippings only
    if islog,
        n = nP*log2(exp(1));
    else
        n = log2(nP);
    end
    
elseif EE && ISE,
    
    % Permutations with sign-flippings
    if islog,
        cte = nP - log(sqrt(pi));
        n   = cte./lambertw(2*cte/exp(1))-.5;
    else
        cte = log(nP/sqrt(pi));
        n   = cte./lambertw(2*cte/exp(1))-.5;
    end
end
