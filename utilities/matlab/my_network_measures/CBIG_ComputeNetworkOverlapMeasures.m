function hubs = CBIG_ComputeNetworkOverlapMeasures(prob_mat, type, threshold)

% hubs = CBIG_ComputeNetworkOverlapMeasures(prob_mat, type, threshold)
%
% prob_mat:  T x N  (probability of location n appears in topic t)
% type:      'global_max', 'entropy', 'vertex_threshold', 'global_thresh',
%            or 'process_cumsum'
% threshold: threshold for types of 'global_max', 'vertex_threshold',
%            'global_thresh', and 'process_cumsum'
%
% 'global_max' gives the entries in prob_mat with higher values than
% threshold times the maximum of prob_mat.
%
% 'entropy' gives the entropy of each cortex location.
%
% 'vertex_threshold' gives the entries in normalized prob_mat (within each
% vertex) with values higher than threshold.
%
% 'global_thresh' gives the entries in prob_mat with higher values than
% threshold.
%
% 'process_cumsum' sorts the vertices descendingly within each topic. The
% vertices with higher cumulative sum of probability are set to 1,
% otherwise 0.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% beta is T x 20484


cortex_label = (sum(prob_mat, 1) ~= 0);
beta = prob_mat(:, cortex_label);
hubs = zeros(1, size(prob_mat, 2));
T = size(beta, 1);
if(strcmp(type, 'global_max'))
    
    beta = (beta > threshold*max(beta(:)));
    hubs(cortex_label) = hubs(cortex_label) + sum(beta, 1);
    
elseif(strcmp(type, 'entropy'))
    
    beta = bsxfun(@times, beta, 1./(squeeze(sum(beta, 1))+eps));
    beta = -sum(beta .* log2(beta), 1); %entropy
    hubs(cortex_label) = hubs(cortex_label) + beta;
    
elseif(strcmp(type, 'vertex_threshold'))
    
    beta = bsxfun(@times, beta, 1./(squeeze(sum(beta, 1))+eps));
    hubs(cortex_label) = hubs(cortex_label) + sum(beta > threshold, 1);
    
elseif(strcmp(type, 'global_thresh'))
    
    beta = (beta > threshold);
    hubs(cortex_label) = hubs(cortex_label) + sum(beta, 1);
    
elseif(strcmp(type, 'process_cumsum'))
    
    beta = beta';
    tmp = sort(beta, 1, 'descend');
    tmp_sum = cumsum(tmp, 1);
    for i = 1:T
        thresh(i) = tmp(find(tmp_sum(:, i) > threshold, 1), i);
        tmp_vol   = squeeze(beta(:, i));
        tmp_vol(tmp_vol < thresh(i))  = 0;
        tmp_vol(tmp_vol >= thresh(i)) = 1;
        beta(:, i) = tmp_vol;
    end
    mean(thresh)
    hubs(cortex_label) = hubs(cortex_label) + sum(beta, 2)';
end