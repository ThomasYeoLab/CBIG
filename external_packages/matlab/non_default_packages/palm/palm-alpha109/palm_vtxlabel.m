function dpvl = palm_vtxlabel(dpv,fac)
% Label a DPV file (vertexwise data).
% 
% Usage:
% dpvl = palm_vtxlabel(dpv,fac)
% 
% - dpv  : Binary data per vertex to be labelled. The non-zero
%          regions receive an unique identifier.
% - fac  : Face indices (see palm_srfread for details).
% - dpvl : Labelled data per vertex.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Feb/2012 (1st version)
% Sep/2013 (this version)
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

% For the labelled data
dpvl = zeros(size(dpv));

% Identify faces that are entirely (3 vertices)
% within the supra-threshold vertices. These are
% the only that will be labelled.
facm = reshape(dpv(fac(:)),size(fac));
faci = all(facm,2);

% The strategy below is as follows: at each iteration of the
% loop, it sorts the labels assigned to the 3 vertices of a
% face. If the last vertex, after sorting is zero, then all 3
% are zero, so none was yet labelled. All three then receive
% the current label "k". If, however, the 3rd is non-zero,
% then all three should stay with that label that was assigned
% in an earlier iteration. But before changing the assignment,
% change also all other vertices that have the same label as
% the 1st and 2nd vertices to have the same as the 3rd.
% Then increment the counter.

% For each face
k = 1; % label indices (to be incremented)
for f = find(faci)',
    
    % Sort the labels assigned for each face
    flab = sort(dpvl(fac(f,:)));
    if flab(3) == 0,
        
        % If the last is zero, then this is the first time all these
        % vertices are seen, so label them with the current "k".
        dpvl(fac(f,:)) = k;

        % Increment the label counter here, as otherwise the labels
        % are reused.
        k = k + 1;
    else
        % Otherwise, this means that the last (3rd) was seen already
        % and is part of some other label. Assign then this higher label
        % to all other vertices in the surface that may have received,
        % in earlier iterations of the loop, indices that are assigned
        % to the 1st or 2nd vertices of the current face.
        if flab(2) ~= 0;
            dpvl(dpvl == flab(2)) = flab(3);
        end
        if flab(1) ~= 0;
            dpvl(dpvl == flab(1)) = flab(3);
        end
        
        % Finally, mark all other vertices of this face with the
        % same lavel as the 3rd one
        dpvl(fac(f,:)) = flab(3);
    end
end
