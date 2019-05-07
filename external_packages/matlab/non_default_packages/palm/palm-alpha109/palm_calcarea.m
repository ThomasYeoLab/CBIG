function A = palm_calcarea(S,isvtx)
% Compute the vertexwise or facewise surface area of a mesh.
% 
% Usage:
% A = palm_calcarea(S,isvtx)
% 
% S     : A struct with fields S.vtx with the vertex coordinates
%         and S.fac with face indices.
% isvtx : A boolean to indicate whether area is to be computed
%         vertexwise (true) or facewise (false).
% A     : Area, either vertexwise or facewise.
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jan/2011 (1st version)
% Sep/2013 (this version, adapted for PALM)
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

% Number of vertices and faces
nV = size(S.vtx,1);
nF = size(S.fac,1);

% Facewise area
facvtx = [S.vtx(S.fac(:,1),:) S.vtx(S.fac(:,2),:) S.vtx(S.fac(:,3),:)];
facvtx0(:,1:6) = facvtx(:,1:6) - [facvtx(:,7:9) facvtx(:,7:9)];
cp  = cross(facvtx0(:,1:3),facvtx0(:,4:6),2);
A   = sqrt(sum(cp.^2,2))./2;

% Vertexwise area
if isvtx,
    dpv = zeros(nV,1);
    A   = A/3;
    for f = 1:nF,
        dpv(S.fac(f,:)) = dpv(S.fac(f,:)) + A(f);
    end
    A = dpv;
end