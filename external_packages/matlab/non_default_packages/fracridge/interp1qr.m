function yi = interp1qr(x,y,xi)
%% Quicker 1D linear interpolation
% Performs 1D linear interpolation of 'xi' points using 'x' and 'y',
% resulting in 'yi', following the formula yi = y1 + (y2-y1)/(x2-x1)*(xi-x1).
% Returns NaN for values of 'xi' out of range of 'x', and when 'xi' is NaN.
%
% 'x'  is column vector [m x 1], monotonically increasing.
% 'y'  is matrix [m x n], corresponding to 'x'.
% 'xi' is column vector [p x 1], in any order.
% 'yi' is matrix [p x n], corresponding to 'xi'.
%
% Copyright (c) 2013 Jose M. Mier
%

%% Full function description
% Quicker 1D linear interpolation: 'interp1qr'
% Performs 1D linear interpolation of 'xi' points using 'x' and 'y',
% resulting in 'yi', following the formula yi = y1 + (y2-y1)/(x2-x1)*(xi-x1).
%
% It has same functionality as built-in MATLAB function 'interp1q' (see
% MATLAB help for details).
%
% It runs at least 3x faster than 'interp1q' and 8x faster than 'interp1',
% and more than 10x faster as m=length(x) increases (see attached performance
% graph).
%
% As with 'interp1q', this function does no input checking. To work properly
% user has to be aware of the following:
%  - 'x'  must be a monotonically increasing column vector.
%  - 'y'  must be a column vector or matrix with m=length(x) rows.
%  - 'xi' must be a column vector.
%
% As with 'interp1q', if 'y' is a matrix, then the interpolation is performed
% for each column of 'y', in which case 'yi' is p=length(xi) by n=size(y,2).
%
% As with 'interp1q', this function returns NaN for any values of 'xi' that
% lie outside the coordinates in 'x', and when 'xi' is NaN.
%
% This function uses the approach given by Loren Shure (The MathWorks) in
% http://blogs.mathworks.com/loren/2008/08/25/piecewise-linear-interpolation/
%  - Uses the function 'histc' to get the 'xi_pos' vector.
%  - Also uses a small trick to rearrange the linear operation, such that
%    yi = y1 + s*(xi-x1), where s = (y2-y1)/(x2-x1), now becomes
%    yi = y1 + t*(y2-y1), where t = (xi-x1)/(x2-x1), which reduces the need
%    for replicating a couple of matrices and the right hand division
%    operation for 't' is simpler than it was for 's' because it takes place
%    only in one dimension (both 'x' and 'xi' are column vectors).
%
% Acknowledgements: Nils Oberg, Blake Landry, Marcelo H. Garcia,
% the University of Illinois (USA), and the University of Cantabria (Spain).
%
% Author:   Jose M. Mier
% Contact:  jmierlo2@illinois.edu
% Date:     August 2013
% Version:  4
%

%% Function begins

% Size of 'x' and 'y'
m = size(x,1);
n = size(y,2);

% For each 'xi', get the position of the 'x' element bounding it on the left [p x 1]
[~,xi_pos] = histc(xi,x);
xi_pos = max(xi_pos,1);     % To avoid index=0 when xi < x(1)
xi_pos = min(xi_pos,m-1);   % To avoid index=m+1 when xi > x(end).

% 't' matrix [p x 1]
dxi = xi-x(xi_pos);
dx = x(xi_pos+1)-x(xi_pos);
t = dxi./dx;

% Get 'yi'
yi = y(xi_pos,:) + t(:,ones(1,n)).*(y(xi_pos+1,:)-y(xi_pos,:));

% Give NaN to the values of 'yi' corresponding to 'xi' out of the range of 'x'
yi(xi<x(1) | xi>x(end),:) = NaN;

end
