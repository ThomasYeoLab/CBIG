function [M,PPH] = palm_vestread(filename)
% Read an FSL "vest" file, e.g. design matrix,
% t-contrasts or f-contrasts.
%
% M = palm_vestread(filename);
%
% filename : File name.
% M        : Matrix.
% PPH      : Peak-to-peak heights.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2012
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

% Read the whole file at once
fid = fopen(filename,'r');
tmp = textscan(fid,'%s');
fclose(fid);

% Get the number of columns
nW = tmp{1}(find(strcmp(tmp{1},'/NumWaves'))+1);
nW = str2double(nW);

% Get the number of rows
nP = tmp{1}(find(strcmp(tmp{1},'/NumPoints') | strcmp(tmp{1},'/NumContrasts'))+1);
nP = str2double(nP);

% Get the peak-to-peak heights
pos = find(strcmp(tmp{1},'/PPheights'));
PPH = tmp{1}(pos+1:pos+nW);
PPH = str2double(PPH)'; % if there is no PPH, this returns NaN

% Reshape to a matrix
M = str2double(tmp{1}(find(strcmp(tmp{1},'/Matrix'))+1:end));
if isempty(M),
    error('File not in the correct format: %s\n',filename);
end
M = reshape(M,[nW,nP])';
