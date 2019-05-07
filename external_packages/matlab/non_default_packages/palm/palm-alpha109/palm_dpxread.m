function varargout = palm_dpxread(filename)
% Read a curvature file (DPV or DPF), in ASCII format.
% This function is much faster than 'dlmread' for large files,
% and works only in Linux and Mac.
%
% [dpx,crd,idx] = dpxread(filename)
% 
% - dpx contains the values for each vertex or face
% - crd contains the vertex coordinates or face indices
% - idx contains vertex or face sequential index
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Feb/2011
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

% Count the number of lines.
% This won't work on Windows but nobody uses it anyway...
[~,result] = system(sprintf('wc -l %s', filename));
nL = str2double(strtok(result,' '));

% Open and read the whole file
fid = fopen(filename,'r');
dpx0 = fscanf(fid,'%f',nL*5);
fclose(fid);

% Reshape from vector to a matrix and get what matters
dpx0 = reshape(dpx0,[5 nL])';
varargout{1} = dpx0(:,5);
varargout{2} = dpx0(:,2:4);
varargout{3} = dpx0(:,1);
