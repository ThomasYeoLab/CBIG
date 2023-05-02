function r = CBIG_ICCW_ICC_1_1(data_vectors)

% function CBIG_ICCW_ICC_1_1(data_vectors)
%
% This code is adapted from the '1-1' formulation for ICC(c) from 
% Arash Salarian, 2008
%
% Reference: McGraw, K. O., Wong, S. P., "Forming Inferences About
% Some Intraclass Correlation Coefficients", Psychological Methods,
% Vol. 1, No. 1, pp. 30-46, 1996
%
% '1-1': The degree of absolute agreement among measurements made on
%         randomly selected objects. It estimates the correlation of any two
%         measurements.
%
% Inputs:
% - data vectors
%   A matrix of observations. Each row is an object of measurement and
%   each column is a judge or measurement.
%
% Outputs:
% - r
%   A scalar. The ICC value.
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% get dimensions
[n, k] = size(data_vectors);
% mean square between subjects
MSR = var(mean(data_vectors, 2)) * k;
% mean square within subjects
MSW = sum(var(data_vectors,0, 2)) / n;

% ICC computation
r = (MSR - MSW) / (MSR + (k-1)*MSW);

end