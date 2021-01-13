function [assign, confidence] = CBIG_IndCBM_wta(surf_labels, vol2surf_fc, topX)

% [assign, confidence] = CBIG_IndCBM_wta(surf_labels, vol2surf_fc, top)
%
% This function use winner-take-all to assign labels to the cerebellum.
% For each cerebellar voxel, this function first identifies the X highest 
% correlated surface vertices. Then, the most frequently assigned network 
% among the X surface vertices is chosen as the "winner" network (each of 
% these X surface vertices is already assigned to a network according to your 
% individual-specific cerebral cortical parcellation). Lastly, we assign the 
% "winner" network to the cerebellar voxel.
% The confidence of the assignment for each cerebellar voxel is measured by
% 1 - N_second_frequent_network/N_most_frequent_network
%
% Input:
%
%     - surf_labels:
%           2N x 1 vector. Parcellation labels of the cereberal cortex.
%           Left hemisphere and right hemisphere are concatenated. 
%
%     - vol2surf_fc:
%           M x N functional connectivity matrix. 
%           M: cerebellar voxels, same order with the cifti template.
%           N: cerebral cortical vertices, lh and rh concatenated.
%           If this matrix is too large and reading the whole matrix takes
%           too much time and memory, you can also save the matrix 
%           'vol2surf_fc' in a mat file and pass in the path. This function
%           will use matfile to read functional connectivity.
%
%     - topX: (double or string)
%           Number of cerebral cortical vertices considered in
%           winner-take-all algorithm. Default is '100'.
%
% Output:
%
%     - assign: 
%           M x 1 vector. Parcellation labels of the cerebellum.
%
%     - confidence: 
%           M x 1 vector. Confidence of the assignment. 
%
% Example:
% [assign, confidence] = CBIG_IndCBM_wta([lh_labels; rh_labels], vol2surf_fc, '100')
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Change 'topX' into num
if(ischar(topX))
    topX = str2num(topX);
end

if(ischar(vol2surf_fc)) % Path of mat file. Use matfile to read the file,
    volmat = matfile(vol2surf_fc);
    [num_voxel, num_vertex] = size(volmat, 'vol2surf_fc');
else
    volmat.vol2surf_fc = vol2surf_fc; % FC matrix directly passed into this function
    [num_voxel, num_vertex] = size(vol2surf_fc);
    clear vol2surf_fc
end
if(num_vertex ~= length(surf_labels))
    error('vol2surf_fc and surf_labels does not match');
end

% Remove medial wall
surf_mask = (surf_labels~=0);
surf_labels = surf_labels(surf_mask);

assign = zeros(num_voxel, 1);
confidence = zeros(num_voxel, 1);

for i = 1: num_voxel
    corr = volmat.vol2surf_fc(i, :); % Read one row each time.
    corr = corr(surf_mask);
    if(sum(isnan(corr)) == length(surf_mask))
        assign(i) = nan;
        confidence(i) = nan;
    else
        corr(isnan(corr)) = -Inf;
        [~, top_idx] = sort(corr, 'descend');
        
        top_idx = top_idx(1: topX);
        top_labels = surf_labels(top_idx); % top X vertices' label info
        
        [freq, unique_labels] = hist(top_labels, unique(top_labels));
        [~, freq_idx] = sort(freq, 'descend');
        
        % Find which network has most vertices (among top X) correlating to
        % the given voxel.
        most_freq_label = unique_labels(freq_idx(1)); 
        assign(i) = most_freq_label;
        
        % Compute confidence for the assignment
        if(length(unique_labels) < 2)
            confidence(i) = 1;
        else
            second_freq_label = unique_labels(freq_idx(2));
            confidence(i) = 1 - sum(top_labels == second_freq_label) / sum(top_labels == most_freq_label);
        end
    end
end

end
