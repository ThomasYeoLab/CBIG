function Y = palm_gammainc(X,A,tail)
% Quick patch until the "gammainc" function in Octave is fixed
% http://savannah.gnu.org/bugs/?47800
% 
% Inputs and outputs are as in "gammainc", except 
% that A must be scalar.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Apr/2016
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2016 Anderson M. Winkler
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

% For most cases, this will be enough:
Y = gammainc(X,A,tail);

% The patch is really only for Octave:
if palm_isoctave,
    idx = Y == 0;
    if any(idx),
        gfun = @(t)exp(-t).*t.^(A-1);
        if strcmpi(tail,'upper'),
            Y(idx) = arrayfun(@(x)quad(gfun,x,Inf),X(idx));
        else
            Y(idx) = arrayfun(@(x)quad(gfun,0,x),X(idx));
        end
        Y(idx) = bsxfun(@rdivide,Y(idx),gamma(A));
    end
end
