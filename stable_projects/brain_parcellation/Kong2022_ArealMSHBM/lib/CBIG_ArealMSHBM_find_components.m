function [ci, sizes]=CBIG_ArealMSHBM_find_components(avg_mesh, ini_labels)

% [ci, sizes]=CBIG_ArealMSHBM_find_components(avg_mesh, ini_labels)
%
% This script will find all distributed components for parcels of a given parcellation
% label <ini_labels> and reorder the parcellation labels of each components. For example,
% <ini_labels> has two parcels 1 and 2, where parcel 1 has 2 components, parcel 2 has 1
% components. Then output <ci> will have 3 ROIs, two of them belong to parcel 1 in the
% original parcellation <ini_labels>.
%
% Input:
%   - avg_mesh:
%     The mesh structures which can be read by CBIG_read_fslr_surface.m or CBIG_ReadNCAvgMesh.m
%     The mesh should be the same as the parcellation space.
%
%   - ini_labels: (Nx1 or 1xN vector)
%     The input parcellation labels. The medial wall should be defined as 0. 
%
% Output:
%   - ci:
%     The updated pacellation labels. Each component will be assigned with a unique label.
%
%   - sizes:
%     The size of each component.
%
% Example:
%   lh_avg_mesh = CBIG_read_fslr_surface('lh','fs_LR_32k','inflated');
%   [lh_ci, lh_sizes]=CBIG_ArealMSHBM_find_components(lh_avg_mesh, lh_labels);    
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

n = size(avg_mesh.vertices, 2);
neighbors=min(size(avg_mesh.vertexNbors));%usually have 6 neighbos
b = zeros(n*neighbors, 1);
c = zeros(n*neighbors, 1);
d = zeros(n*neighbors, 1);
for i = 1:12 % vertices with fewer neighbors (5 neighbors)
    b((i-1)*neighbors+1:i*neighbors-1)=i;
    c((i-1)*neighbors+1:i*neighbors-1)=avg_mesh.vertexNbors(1:neighbors-1, i);
    d((i-1)*neighbors+1:i*neighbors-1)=(ini_labels(avg_mesh.vertexNbors(1:neighbors-1, i))==ini_labels(i));
end
for i = 13:n % vertices with regual neighbors (6 neighbors)
    b((i-1)*neighbors+1:i*neighbors)=i;
    c((i-1)*neighbors+1:i*neighbors)=avg_mesh.vertexNbors(:, i);
    d((i-1)*neighbors+1:i*neighbors)=(ini_labels(avg_mesh.vertexNbors(:, i))==ini_labels(i));
end
for i = 12:-1:1% account for zero vertices
    b(i*neighbors)=[];
    c(i*neighbors)=[];
    d(i*neighbors)=[];
end
g = sparse(b,c,d);
[ci, sizes] = components(g);

end