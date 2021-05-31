function [lh_labels,rh_labels]=CBIG_ArealMSHBM_BuildTwoVertThickBoundary(lh_avg_mesh, ...
    rh_avg_mesh, lh_orig_labels,rh_orig_labels)

% [lh_labels,rh_labels]=CBIG_ArealMSHBM_BuildTwoVertThickBoundary(lh_avg_mesh, ...
% rh_avg_mesh, lh_orig_labels,rh_orig_labels)
%
% This script will create boundaries for a given parcellation <lh_orig_labels> 
% <rh_orig_labels>. The thickness of boundaries are two vertices. Specifically, vertices
% whose neighbors are assigned with different labels will be considered as boundary
% vertices. The boundary vertices will be denoted as 0 in the output <lh_labels>,
% <rh_labels>. 
%
% Input:
%   - lh_avg_mesh, rh_avg_mesh:
%     The mesh structures which can be read by CBIG_read_fslr_surface.m or CBIG_ReadNCAvgMesh.m
%     The mesh should be the same as the parcellation space.
%
%   - lh_orig_labels, rh_orig_labels: (Nx1 vector)
%     The input parcellation labels. The medial wall should be defined as 0.     
%
% Output:
%   - lh_labels, rh_labels:
%     The output parcellation with boundary vertices. The boundary vertices will be denoted as 0s.
%
% Example:
%   mesh = 'fsaverage6';
%   lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','aparc.annot');
%   rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated','aparc.annot');
%   lh_group_labels = CBIG_read_annotation('lh.Schaefer2018_400Parcels_17Networks_order.annot']);
%   rh_group_labels = CBIG_read_annotation('rh.Schaefer2018_400Parcels_17Networks_order.annot']);
%   lh_orig_labels = lh_group_labels - 1;
%   rh_orig_labels = rh_group_labels - 1;
%   [lh_labels,rh_labels]=CBIG_ArealMSHBM_BuildTwoVertThickBoundary(lh_avg_mesh, ...
%       rh_avg_mesh, lh_orig_labels,rh_orig_labels);
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


lh_neigh_lables=lh_avg_mesh.vertexNbors;
rh_neigh_lables=rh_avg_mesh.vertexNbors;
lh_neigh_lables(lh_neigh_lables~=0)=lh_orig_labels(lh_neigh_lables(lh_neigh_lables~=0));
rh_neigh_lables(rh_neigh_lables~=0)=rh_orig_labels(rh_neigh_lables(rh_neigh_lables~=0));

lh_full_labels=[lh_orig_labels';lh_neigh_lables];
rh_full_labels=[rh_orig_labels';rh_neigh_lables];
lh_temp=bsxfun(@minus,lh_full_labels,lh_full_labels(1,:));
rh_temp=bsxfun(@minus,rh_full_labels,rh_full_labels(1,:));
lh_temp(lh_full_labels==0)=0;
rh_temp(rh_full_labels==0)=0;
lh_labels=lh_orig_labels;
rh_labels=rh_orig_labels;
lh_labels(sum(lh_temp)~=0)=0;
rh_labels(sum(rh_temp)~=0)=0;