function icc = CBIG_fisher_icc(data1, data2)

% icc = CBIG_fisher_icc(data1, data2)
% Compute fisher's intraclass correlation between two vectors.(This function will exclude NAN)
% 
%   cc = CBIG_fisher_icc(data1, data2)
%   Input:
%       data1: N x 1 vector
%       data2: N x 1 vector
%   Output:
%       icc  : 1 x 1 scalar
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(data1,2) ~= 1
    error('Input argument ''lh_data'' should be a column vector');
end
if size(data2,2) ~= 1
    error('Input argument ''rh_data'' should be a column vector');
end

data = [data1; data2];
data_mean = nanmean(data);
data_var  = nanvar(data, 1);

icc = sum((data1 - data_mean).*(data2 - data_mean)) /(numel(data1) * data_var);  

