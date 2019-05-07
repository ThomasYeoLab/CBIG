function [vtx,fac] = palm_srfread(filename)
% Read a surface file, in ASCII format.
% 
% [vtx,fac] = palm_srfread('filename');
% 
% - vtx contains the coordinates (x,y,z), one vertex per row
% - fac contains the indices for the three vertices of each face
% 
% The indices for the vertices start at 1, not 0.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jan/2011
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

fid = fopen(filename,'r');
fgets(fid); % ignore 1st line
nV  = fscanf(fid,'%d',1); % num of vertices
nF  = fscanf(fid,'%d',1); % num of faces
vtx = fscanf(fid,'%f',nV*4);
fac = fscanf(fid,'%d',nF*4);
fclose(fid);
vtx = reshape(vtx,[4 nV])';
vtx(:,4) = [];
fac = reshape(fac,[4 nF])';
fac(:,4) = [];
fac = fac + 1; % indices start at 1, not 0 in OCTAVE/MATLAB.
