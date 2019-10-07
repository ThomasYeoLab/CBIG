function CBIG_ASDf_plotKmeans(C, file_name, scalelim)
% CBIG_ASDf_plotKmeans(C, file_name, scalelim)
% 
% Plot k-means cluster centroids.
%
% Input:
%     - C:
%           Cluster centroids obtained from k-means algorithm, number of
%           clusters x number of unique ROI-ROI pairs
%     - file_name:
%           File name prefix to save the plot.
%     - scalelim:
%           Scale range for the plot. E.g., [-1.8e-3 1.8e-3].
%
% Example:
%       CBIG_ASDf_plotKmeans(C, '~/output/kmeansCluster', [-1.6e-3 1.6e-3])
% 
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_clusters = size(C,1);
num_ROIs = 419;

if size(C,2) ~= (num_ROIs * (num_ROIs - 1) / 2)
    error('Only configured for 419x419 ROI.');
end

for idx = 1:num_clusters
    centroid = C(idx,:);
    save([file_name '_cluster' num2str(idx) '.mat'], 'centroid');

    % Get back to 419x419 matrix
    index = 0;
    kmeans_c = zeros(num_ROIs, num_ROIs);
    for j = 1:(num_ROIs-1)
        for i = (j+1):num_ROIs
            index = index + 1;
            kmeans_c(i,j) = centroid(index);
        end
        kmeans_c(i,i) = 0; % Diagonal set to 0
    end
    for i = 1:num_ROIs
        for j = (i+1):num_ROIs
            kmeans_c(i,j) = kmeans_c(j,i);
        end
    end
    
    % Plot the matrix
    if nargin <= 2 || isempty(scalelim)
        scalelim = [-max(abs(kmeans_c(:))) max(abs(kmeans_c(:)))];
    end
    CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(kmeans_c, scalelim, [file_name '_cluster' num2str(idx)]);
  
end