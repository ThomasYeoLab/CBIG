function y = palm_isoctave
% Test whether the function is running in OCTAVE or not.
% Retuns 1 (true) if yes, 0 (false) otherwise.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2012
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

persistent palm_isoct;
if isempty(palm_isoct),
    palm_isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end
y = palm_isoct;
