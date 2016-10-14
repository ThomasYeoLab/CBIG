function [m,ssq] = issq(x,n,m0,ssq0)
% Increment the mean and the sum of squares after a new
% observation is included in the sample.
% To obtain the variance, divide the ssq by n or by n-1.
%
% Usage:
% [m,ssq] = issq(x,n,m0,ssq0)
%
% Inputs:
% x    = New observation
% n    = Observation sequential number (n > 1)
% m0   = Previous mean
% ssq0 = Previous sum of squares
%
% Outputs:
% m    = Updated mean
% ssq  = Updated sum of squares
%
% Reference:
% * Welford, B.P. Note on a method for calculating corrected sums of
%   squares and products. Technometrics (1962) 4(3):419-20
% 
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Mar/2011

% Accept arguments
if nargin < 3,
    error('Insufficient arguments.')

elseif n <= 1 || round(n) ~= n,
    error('N has to be integer and > 1.')

end

% Compute the incremented mean
m = ((n-1)*m0 + x)/n;

% Compute the incremented sum of squares
if nargout == 2,
    if n == 1,
        ssq = zeros(size(x));

    elseif nargin == 4,
        ssq = ssq0 + ((n-1)*(x-m0).^2)/n;

    else
        error('Previous SSQ not supplied.');

    end
end
