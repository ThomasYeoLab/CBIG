function stability = CBIG_ComputeClusterStability(labels)

% stability = CBIG_ComputeClusterStability(labels)
% 
% Compute how stable the parcellation labels are in multiple trials.
% Stability is defined by Dice overlap.
%
% assume labels are num_vertices x # trials
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



N = size(labels, 1);
T = size(labels, 2);
K = max(labels(:, 1));
stability = zeros(N, 1);

tic
dice_mat = zeros(K, K, T);
for t = 1:T % foreach reference
       
    disp(num2str(t)); 
    
    ref_labels = labels(:, t);
    for i = 1:K
        ref_vec = double(ref_labels == i);  
        for j = 1:K
            neighbors = double(labels == j);
            innerprod = sum(bsxfun(@times, neighbors, ref_vec), 1);
            dice_mat(i, j, :) = 2*innerprod./(sum(neighbors, 1) + sum(ref_vec));
        end
    end
    
    [a, c] = meshgrid(labels(:, t), 1:T);
    a = a(:);
    c = c(:);
    ind_vec = sub2ind(size(dice_mat), a, reshape(labels', numel(labels), 1), c);
    dice = reshape(dice_mat(ind_vec), T, N);
    stability = stability + mean(dice, 1)';
end
stability = stability/T;
toc
    
    
    
    
