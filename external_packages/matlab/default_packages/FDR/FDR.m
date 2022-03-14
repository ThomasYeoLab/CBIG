function [ ind, thres ] = FDR(p_list, alpha)

% This code has been modified from the original version written by
% Edden Gerber, lab of Leon Y. Deouell, 2012. This code computes the 
% False Discovery Rate according to Benjamini and Hochberg (1995). 
% 
% Inputs: 
% - p_list 
%   A vector. This should be the complete list of p values calculated 
%   from your experiment.
%
% - alpha 
%   The desired alpha threshold for FDR correction. Default: 0.05
%
%
% Outputs:
% - ind
%   The indexes of significant p-values within p_list.
%
% - thres 
%   The p-value which served as the actual threshold in this test. 
% 
% Adapted by Leon Ooi and CBIG.

% find number of p-values for correction
n_vals = length(p_list);
num_tests = n_vals; 

% set alpha to 0.05 if not passed in
if nargin < 2
    alpha = 0.05;
end

% calculate corrected value
p_sorted = sort(p_list,'descend');
comp = (num_tests:-1:1)/num_tests * alpha;
comp = comp((end-n_vals+1):end);

% return values that pass FDR
i = find(p_sorted <= comp,1,'first');
if isempty(i)
    thres = 0;
else
    thres = p_sorted(i);
end
ind = find(p_list<=thres);

end

