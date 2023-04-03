function [ci, sizes] = CBIG_hMRF_generate_components_one_hemi(hemi_avg_mesh, hemi_full_label)
% [ci, sizes] = CBIG_hMRF_generate_components_one_hemi(hemi_avg_mesh, hemi_full_label)
%
% This function computes connected components for a labeled vector.

% For the notations below:
% N = no of vertices per hemisphere; 

% Input
%   - hemi_avg_mesh: (struct)
%     The struct containing mesh information for a single hemisphere.
%   - hemi_full_label: (matrix)
%     Nx1 parcellation label for a single hemisphere.
% Output
%   - ci (matrix)
%     A Nx1 vector indicating the segmented connected components. In short, each connected component will share
%     one unique label.
%   - sizes (matrix)
%     A Nx1 vector indicating the sizes of the segmented connected components.
%
% Example
%   - [ci, sizes] = CBIG_hMRF_generate_components_one_hemi(hemi_avg_mesh, hemi_full_label)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ini_labels = hemi_full_label;
n = size(hemi_avg_mesh.vertices, 2);
neighbors = min(size(hemi_avg_mesh.vertexNbors)); % usually 6
b = zeros(n*neighbors, 1);
c = zeros(n*neighbors, 1);
d = zeros(n*neighbors, 1);
for i = 1: 12 % vertices with fewer neighbors
    b((i-1)*neighbors+1: i*neighbors-1) = i;
    c((i-1)*neighbors+1: i*neighbors-1) = hemi_avg_mesh.vertexNbors(1:neighbors-1, i);
    d((i-1)*neighbors+1: i*neighbors-1) = ...
        (ini_labels(hemi_avg_mesh.vertexNbors(1:neighbors-1, i)) == ini_labels(i));
end
for i = 13: n % vertices with regular neighbors
    b((i-1)*neighbors+1: i*neighbors) = i;
    c((i-1)*neighbors+1: i*neighbors) = hemi_avg_mesh.vertexNbors(:, i);
    d((i-1)*neighbors+1: i*neighbors) = (ini_labels(hemi_avg_mesh.vertexNbors(:, i)) == ini_labels(i));
end
for i = 12: -1: 1% account for zero vertices
    b(i*neighbors) = [];
    c(i*neighbors) = [];
    d(i*neighbors) = [];
end
g = sparse(b, c, d);
[ci, sizes] = components(g);
end
