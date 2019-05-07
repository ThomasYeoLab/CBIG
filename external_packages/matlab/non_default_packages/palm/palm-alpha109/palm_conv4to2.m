function Y2d = palm_conv4to2(Y4d)
% Convert a 4D dataset (x,y,z,t) into a 2D array (t,x*y*z)
% that can be used in a multivariate GLM
% 
% Usage:
% Y2d = conv4to2(Y4d);
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

% tmp = permute(Y4d,[4 1 2 3]);
% siz = size(tmp);
% Y2d = reshape(tmp,[size(tmp,1) prod(siz(2:end))]);

Y2d = reshape(Y4d,numel(Y4d)/size(Y4d,4),size(Y4d,4))';

