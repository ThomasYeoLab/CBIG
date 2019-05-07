function B = palm_incrbin(B)
% Increment a binary number by 1. 
% 
% Usage:
% Bi = palm_incrbin(B)
%
% B  : A logical vector, with the least
%      significant digits first.
% Bi : Incremented B.
% 
% If the highest has been reached, Bi = B.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014
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

k = find(~B,1);
B(1:k) = ~B(1:k);
