function palm_msetwrite(varargin)
% Write a .mset file, i.e. an ASCII file
% containing multiple matrices (2D arrays).
%
% palm_msetwrite(filename,Mset)
% palm_msetwrite(filename,M1,M2,M3,...)
%
% filename : File name to be created.
% Mset     : Cell array with the arrays.
% M1, ...  : Contrasts as individual arrays.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2015
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

% Some basic argument checking
if nargin < 2,
    error('At least 2 input arguments must be supplied');
elseif ~ischar(varargin{1}),
    error('The first argument must be a filename');
end

% Set of contrasts to be saved
if nargin == 2 && iscell(varargin{2}),
    Mset = varargin{2};
else
    Mset = varargin(2:end);
end
nM = numel(Mset);

% Before saving, make sure all are 2D arrays
for m = 1:nM,
    nd = ndims(Mset{m});
    if nd > 2,
        error('Inputs must be 2D arrays (ndims of argument #%d is %d',m+1,nd);
    end
end

% Write to the disk
fid = fopen(varargin{1},'w');
for m = 1:nM,
    fstr = horzcat('%0.6e',repmat('\t%0.6e',1,size(Mset{m},2)-1),'\n');
    fprintf(fid,'Matrix %d %d\n',size(Mset{m}));
    fprintf(fid,fstr,Mset{m}');
    fprintf(fid,'\n');
end
fclose(fid);
