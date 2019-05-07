function palm_ptree2dot(Ptree,dotfile)
% Create a DOT file from a permutation tree, that can be used
% with GraphViz for visualisation.
%
% Usage:
% palm_ptree2dot(Ptree,dotfile)
%
% Ptree     : Permutation (dependence) tree.
% dotfile   : DOT file to be created.
%
% Reference:
% * Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM.
%   Multi-level block permutation. Neuroimage. 2015;123:253-68.
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

graphname = 'Ptree';
ntop = 'x1';
fid = fopen(dotfile,'w');
fprintf(fid,'graph %s {\n',graphname);
if size(Ptree,2) > 1,
    if isnan(Ptree{1,1}(1)),
        nclr = 'red';
    else
        nclr = 'blue';
    end
    fprintf(fid,'%s [label="" shape=point color=%s fixedsize=true width=1 height=1 fontsize=20 penwidth=3];\n',ntop,nclr);
    downstream(fid,ntop,Ptree{1,3});
end
fprintf(fid,'}');
fclose(fid);

% ==============================================================
function downstream(fid,ntop,Ptree)
for u = 1:size(Ptree,1),
    
    % Current node
    ncur    = sprintf('%sx%s',ntop,num2str(u));
    
    % Print node
    inancur = isnan(Ptree{u,1}(1));
    if size(Ptree,2) == 1,
        nlab = sprintf('"%d"',Ptree{u,1});
        nshp = 'circle';
        nclr = 'black';
        nwid = '1';
        nhei = '1';
    elseif size(Ptree{u,3},1) == 1,
        nlab = '""';
        nshp = 'point';
        nclr = 'black';
        nwid = '.2';
        nhei = '.2';
    else
        nlab = '""';
        nshp = 'point';
        if inancur,
            nclr = 'red';
        else
            nclr = 'blue';
        end
        nwid = '1';
        nhei = '1';
    end
    fprintf(fid,'%s [label=%s shape=%s color=%s fixedsize=true width=%s height=%s fontsize=20 penwidth=3];\n',...
        ncur,nlab,nshp,nclr,nwid,nhei);
    
    % Print edge
    fprintf(fid,'%s -- %s [color=black penwidth=6];\n',ntop,ncur);
    
    % Keep going down more levels as needed
    if size(Ptree,2) > 1,
        downstream(fid,ncur,Ptree{u,3});
    end
end
