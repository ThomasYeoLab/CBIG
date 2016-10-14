function [low_bound, upp_bound] = CBIG_bootstrap_confidence_interval(input_data, alpha, nsamples, stat_functor, varargin)

% Calculate the confidence interval from bootstrap distribution.
% 
%   [low_bound, upp_bound] = CBIG_bootstrap_confidence_interval(input_data, alpha, nsamples, stat_functor, varargin)
%   Input:
%       input_data  : some dimensional matrix, where last dimension is the number of data points.
%       alpha       : significance level
%       nsamples    : number of resamples with replacement  
%       stat_functor: function name, input_data and varargin are the input of the function
%       varargin    : input arguments of the above function
%   Output:
%       low_bound   : lower bound of confidence interval
%       upp_bound   : upper bound of confidence interval
% 
%   Example:
%       [low_bound, upp_bound] = CBIG_bootstrap_confidence_interval([1 2 3; 4 5 6], 0.05, 1000, @prod);
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% set the state of the random seed. Each time you run this code, you will get the same random number.
rand('state', 5489); 

% get the length of last dimension for input_data
size_vec = size(input_data);
N = size(input_data, length(size_vec)); 

% rearrange to two dimensional matrix
input_data = reshape(input_data, prod(size_vec(1:length(size_vec)-1)), N);

% bootstrapping!
for i = 1:nsamples
    data = reshape(input_data(:, randsample(N, N, 1)), size_vec);
    tmp = feval(stat_functor, data, varargin{:});
    
    if(i == 1)
       stat = zeros(nsamples, length(tmp)); 
    else
       stat(i, :) = tmp; 
    end
end

% compute confidence intervals
stat = sortrows(stat);
low_bound = stat(round(alpha/2 * nsamples), :);
upp_bound = stat(round(nsamples - alpha/2 * nsamples), :);
