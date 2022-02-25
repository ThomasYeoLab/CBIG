function [ ind, thres ] = FDR( p_list, alpha, corrected )
% Computes the False Discovery Rate according to Benjamini and Hochberg (1995). 
% 
% Inputs: 
% p_list - list of p values
% alpha - the desired alpha threshold. Default: 0.05
% corrected - set to true if correction for dependencies is to be applied, according to Benjamini
% and Yekutieli (2001) (this is probably not the common case). 
%
% outputs:
% ind - the indexes of significant p-values within p_list
% thres - the p-value which served as the actual threshold in this test. 
% 
% Written by Edden Gerber, lab of Leon Y. Deouell, 2012
% Please send bug reports and requsts to edden.gerber@gmail.com
%
n_vals = length(p_list);
num_tests = n_vals; % there was some reason that in some cases you may want to set this to
% a lower value, but I don't remember what it was. 

if nargin < 2
    alpha = 0.05;
end

if nargin < 3
    corrected = false;
end

p_sorted = sort(p_list,'descend');

if corrected
    comp = (num_tests:-1:1)/num_tests * alpha / sum((1:num_tests)/num_tests);
else
    comp = (num_tests:-1:1)/num_tests * alpha;
end


comp = comp((end-n_vals+1):end);

i = find(p_sorted <= comp,1,'first');

if isempty(i)
    thres = 0;
else
    thres = p_sorted(i);
end

ind = find(p_list<=thres);

end

