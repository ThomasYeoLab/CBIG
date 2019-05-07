function dpfl = palm_faclabel(dpf,fac)
% Label a DPF file (facewise data).
% 
% Usage:
% dpfl = palm_faclabel(dpf,fac)
% 
% - dpf  : Binary data per face to be labelled. The non-zero
%          regions receive an unique identifier.
% - fac  : Face indices (see palm_srfread for details).
% - dpfl : Labelled data per vertex.
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

% Ignore faces that are not in the mask
dpfl = zeros(size(dpf));
facm = bsxfun(@times,fac,dpf);

% For each face
k = 1; % label indices (to be incremented)
for f = find(dpf)',
    
    % Find other faces that share vertices with
    % the current one
    neifacidx = ...
        facm == fac(f,1) | ...
        facm == fac(f,2) | ...
        facm == fac(f,3);
    
    % Only 2+ shared vertices, not 1 (i.e., edge in common)
    neifacidx = sum(neifacidx,2) >= 2;
    
    % Number of neighbours
    numneigh  = sum(neifacidx);
    
    % If this isn't an isolated face
    if numneigh > 1,
        
        % Sort the labels assigned to the neighbours
        flab = sort(dpfl(neifacidx));
        if flab(numneigh) == 0,
            
            % If the last (highest) label of the neighbours is
            % zero, this means neither neighbour has received a
            % label (yet, perhaps never). So use a new label, 
            % and increment the counter.
            dpfl(neifacidx) = k;
            k = k + 1;
        else
            
            % If one or more neighbours already have a label, assign
            % the highest to all the others, including all other
            % faces across the surface that may have received those.
            for n = numneigh:-1:1,
                if flab(n) ~= 0,
                    dpfl(dpfl == flab(n)) = flab(numneigh);
                end
            end
            
            % Only then assign the highest label to the current face
            dpfl(f) = flab(numneigh);
        end
        
    else
        
        % If this is an isolated face, label it and increment the counter.
        dpfl(f) = k;
        k = k + 1;
    end
end
