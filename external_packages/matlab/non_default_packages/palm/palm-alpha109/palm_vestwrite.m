function palm_vestwrite(filename,M)
% Write an FSL "vest" file, e.g. design matrix,
% t-contrasts or f-contrasts.
%
% palm_vestwrite(filename,M);
%
% filename : File name.
% M        : Matrix.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
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

% Number of waves, points and peak-to-peak heights
[nP,nW] = size(M);
PPH = max(M,[],1) - min(M,[],1);

% Formatting string
fstr = horzcat('%0.6e',repmat('\t%0.6e',1,nW-1),'\n');

% Write to the disk
fid = fopen(filename,'w');
fprintf(fid,'/NumWaves\t%d\n',nW);
fprintf(fid,'/NumPoints\t%d\n',nP);
fprintf(fid,horzcat('/PPHeights\t',fstr),PPH);
fprintf(fid,'\n/Matrix\n');
fprintf(fid,fstr,M');
fclose(fid);
