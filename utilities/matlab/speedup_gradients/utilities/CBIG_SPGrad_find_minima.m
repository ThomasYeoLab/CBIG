function minimametric = CBIG_SPGrad_find_minima(data,K_neighbors)

% minimametric = CBIG_gMSHBM_find_minima(data,K_neighbors)
%
% This script find the local minima of the smoothed gradient map
%
% Input:
%     - data:
%       the smoothed gradient map.
%
%     - K_neighbors:
%       the K nearest neighbors for a mesh.
%
% Output:
%     - minimametric:
%       the local minima of the smoothed gradient map
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%fake a datapoint for NaN neighbors
K_neighbors(:,1) = [];
K_neighbors(isnan(K_neighbors)) = size(K_neighbors,1) + 1;
K_neighbors = [K_neighbors; ones(1,size(K_neighbors,2))*(size(K_neighbors,1) + 1)];
data = [data; Inf(1,size(data,2))];

for i = 1:size(K_neighbors,2)
    tmp = (data - data(K_neighbors(:,i),:)) < 0;
    if(i ==1)
        minimametric = tmp;
    else
        minimametric = tmp & minimametric;
    end
end
minimametric(end,:) = [];    
