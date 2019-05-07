function Br = palm_reindex(varargin)
% Reindexes blocks so that each is numbered natural numbers,
% starting from 1. There are two pure possibilities: the numbering
% goes continuously for each level, crossing block definitions,
% or it restarts at 1 for each new block. A third method is a
% mixture of both, i.e., it restarts at each block, but at the
% last block it is continuous, crossing blocks.
% Finally (and the default), it is also possible to add one
% extra column to treat the simplification in which the last
% level is omited for whole-block permutation ('fixleaves').
%
% Usage:
% Br = palm_reindex(B,method)
% 
% B     : Block definitions (multilevel).
% meth  : Method for reindexing: 'continuous', 'restart',
%         'mixed' or 'fixleaves' as described above.
%         Default: 'fixleaves'.
% Br    : Reindexed block definitions.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Oct/2013
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

% Default
meth = 'fixleaves';

% Take care of input and output vars
B  = varargin{1};
Br = zeros(size(B));
if nargin > 1,
    meth = varargin{2};
end

switch lower(meth),
    
    case 'continuous',
        
        % Renumber using sequential natural numbers starting
        % from 1 and crossing blocks. The first level is
        % treated differently as it doesn't have a
        % "previous" block.
        U = unique(B(:,1));
        for u = 1:numel(U),
            idx = B(:,1) == U(u);
            Br(idx,1) = u*sign(U(u));
        end
        for b = 2:size(B,2), % 2nd onwards
            Bb = B(:,b);
            Bp = Br(:,b-1); % previous
            Up = unique(Bp);
            cnt = 1;
            for up = 1:numel(Up),
                idxp = Bp == Up(up);
                U = unique(Bb(idxp));
                for u = 1:numel(U),
                    idx = (Bb == U(u)) & idxp;
                    Br(idx,b) = cnt*sign(U(u));
                    cnt = cnt + 1;
                end
            end
        end
        
    case 'restart',
        
        % Renumber using sequential natural numbers
        % starting from 1 but never crossing blocks,
        % restarting instead for each block.
        Br = renumber(B);
        
    case 'mixed'
        
        % This mixes both above
        Ba = palm_reindex(B,'restart');
        Bb = palm_reindex(B,'continuous');
        Br = horzcat(Ba(:,1:end-1),Bb(:,end));
        
    case 'fixleaves',
        
        % Renumber using sequential natural numbers
        % starting from 1 but never crossing blocks,
        % restarting instead for each block.
        [Br,addcol] = renumber(B);
        if addcol,
            Br = horzcat(Br,(1:size(Br,1))');
            Br = renumber(Br);
        end
        
    otherwise
        error('Unknown method: %s',meth);
end

% ==============================================
function [Br,addcol] = renumber(B)
% Note that this runs recursively.

B1 = B(:,1);
U = unique(B1);
addcolvec = false(size(U));
nU = numel(U);
Br = zeros(size(B));
for u = 1:nU,
    idx = B1 == U(u);
    Br(idx,1) = u*sign(U(u));
    if size(B,2) > 1,
        [Br(idx,2:end),addcolvec(u)] = renumber(B(idx,2:end));
    elseif sum(idx) > 1,
        addcol = true;
        Br(idx) = -abs(B(idx));
    else
        addcol = false;
    end
end

if size(B,2) > 1,
    addcol = any(addcolvec);
end
