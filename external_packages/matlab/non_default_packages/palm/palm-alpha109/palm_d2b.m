function b = palm_d2b(d,n)
% Converts decimal integers to binary.
% 
% b = palm_d2b(d,n)
%
% d : Decimal value(s).
% n : Word size (number of columns of b).
% b : Binary value(s).
%     The least significant go last.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2013
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

% Although this function is tiny, it's used more than
% once, so it's better to have as a separate file.

[~,e] = log2(max(d));
b = rem(floor(d(:)*pow2(1-max(n,e):0)),2);