function CBIG_pFIC_draw_surface_maps(input_vertex, fig_title, fig_save_path, varargin)

% CBIG_pFIC_draw_surface_maps(input_vertex, fig_title)
% This funciton plots a surface given a vector of numbers in fsaverage5 space at the
% vertex-level (exclude the medial-wall)
% Input:
%   - input_vertex: a 18715-by-1 vector. fsaverage5 space has 20484 vertices
%   in total, but 1769 of them are medial walls. Excluding them leaves
%   18715 vertices.
%   - fig_title: title of the figure, note that this will also be used as
%   the file name of an intermediate output, so do not leave it empty
%   - fig_save_path: path to save out the figures
%   - varargin: empty or 'cognition', if varargin == 'cognition', flip the
%   colormap. this is to be consistent with the age effect results in that cyan
%   indicates a larger effect (steeper slope or higher cohen d) and magenta
%   indicates a smaller effect (gentler slope or lower cohen d). (See CBIG_pFIC_plot_figure_5.m Figure 5 (D) & (F))
% Example:
% CBIG_pFIC_draw_surface_maps(vertex_level_EI_ratio_in_fsaveage5, 'EI_ratio', '/home/')
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% create a temporary working directory
tmp_dir = './tmp/';
if ~exist(tmp_dir ,'dir')
    mkdir(tmp_dir);
end

%% upsample the labels to fsaverage6 space then project it to fsLR 164k space
lh_avg_mesh_fs5 = CBIG_ReadNCAvgMesh('lh','fsaverage5','inflated', 'cortex');
rh_avg_mesh_fs5 = CBIG_ReadNCAvgMesh('rh','fsaverage5','inflated', 'cortex');

lh_avg_mesh_fs6 = CBIG_ReadNCAvgMesh('lh','fsaverage6','inflated', 'cortex');
rh_avg_mesh_fs6 = CBIG_ReadNCAvgMesh('rh','fsaverage6','inflated', 'cortex');

lh_mask = (lh_avg_mesh_fs5.MARS_label == 2); % 1 for medial wall, 2 for cortex
rh_mask = (rh_avg_mesh_fs5.MARS_label == 2); 

% reorder labels to continuous numbers
isNeg = all(input_vertex <= 0);
unique_vertex = unique(input_vertex);
for i = 1:length(unique_vertex)
    input_vertex(input_vertex == unique_vertex(i)) = i;
end

lh_label = input_vertex(1:sum(lh_mask));
rh_label = input_vertex(sum(lh_mask)+1:end);

lh_full_label_fs5 = zeros(size(lh_mask));
lh_full_label_fs5(lh_mask) = lh_label;

rh_full_label_fs5 = zeros(size(rh_mask));
rh_full_label_fs5(rh_mask) = rh_label;

% delete existing .mat file
if exist(fullfile(tmp_dir, [fig_title '.mat']),'file')
    delete(fullfile(tmp_dir, [fig_title '.mat']));  
end

% upsample labels to fsaverage6
lh_full_label_fs6 = MARS_NNInterpolate_kdTree(lh_avg_mesh_fs6.vertices, lh_avg_mesh_fs5, lh_full_label_fs5);
rh_full_label_fs6 = MARS_NNInterpolate_kdTree(rh_avg_mesh_fs6.vertices, rh_avg_mesh_fs5, rh_full_label_fs5);
% project to fsLR164k
[~, ~, lh_label_fsLR_164k, rh_label_fsLR_164k] = CBIG_project_fsaverage2fsLR(lh_full_label_fs6,...
    rh_full_label_fs6, 'fsaverage6', 'label', ['~/tmp_' num2str(rand(1))]);
save(fullfile(tmp_dir, [fig_title '.mat']), 'lh_label_fsLR_164k', 'rh_label_fsLR_164k');

lh_annot = 'fsLR_164k_lh_annot.annot';
rh_annot = 'fsLR_164k_rh_annot.annot';

ct = colormap(cool(length(unique(input_vertex))));
ct = [0,0,0; ct]; % in annot files, we usually use 0 to represent the medial wall
if ~isNeg
    ct(2, :) = [0 0 0];
else 
    ct(end, :) = [0 0 0];
end
ct_large = round(ct*255);

% the masked ROIs are plotted as silver grey. This is to differentiate from
% the medial wall, which is plotted in black
% (for example, in the case of Alprazolam E/I ratio contrsat, temporal pole (ROI 34) does not show a significant E/I 
% difference after FDR. Thus, the temporal pole is masked out when plotting the surface map and it is shown in 
% silver grey)
if sum(isinf(unique_vertex)) ~= 0
    ct_large(end, :) = [192 192 192]; 
end

% when plottiing coginition effect cohen d results, flip the color bar
% this is to be consistent with the age effect results in that cyan
% indicates a larger effect (steeper slope or higher cohen d) and magenta
% indicates a smaller effect (gentler slope or lower cohen d).
if strcmp(varargin, 'cognition')
   ct_large(3:end, :) = flip(ct_large(3:end, :)); 
end
load('yeo7network_fsLR164k_label.mat', 'lh_yeo_label_fsLR_164k', 'rh_yeo_label_fsLR_164k');

mesh_name = 'fs_LR_164k';
lh_avg_mesh_fs164k = CBIG_read_fslr_surface('lh', mesh_name, 'very_inflated');
rh_avg_mesh_fs164k = CBIG_read_fslr_surface('rh', mesh_name, 'very_inflated');

[yeo_lh_boundary, yeo_rh_boundary] = ...
    CBIG_pFIC_BuildTwoVertThickBoundary(lh_avg_mesh_fs164k, rh_avg_mesh_fs164k, ...
    lh_yeo_label_fsLR_164k, rh_yeo_label_fsLR_164k);
lh_label_fsLR_164k(yeo_lh_boundary == 0) = 0;
rh_label_fsLR_164k(yeo_rh_boundary == 0) = 0;

CBIG_pFIC_parcellation_to_annot(lh_label_fsLR_164k, rh_label_fsLR_164k, lh_annot, rh_annot, ct_large);

%% plot with freesurfer
command = ['freeview -f $CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/lh.very_inflated:annot=' ...
    lh_annot ':color=5,25,5:edgethickness=0 --viewport 3d -cam dolly 1.6 -ss ' fig_title '_lh_lateral_tmp.png'];
system(command);
command = ['freeview -f $CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/lh.very_inflated:annot=' ...
    lh_annot ':color=5,25,5:edgethickness=0 --viewport 3d -cam dolly 1.6 Azimuth 180 -ss ' ...
    fig_title '_lh_medial_tmp.png'];
system(command);
command = ['freeview -f $CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/rh.very_inflated:annot=' ...
    rh_annot ':color=5,25,5:edgethickness=0 --viewport 3d -cam dolly 1.6 -ss ' fig_title '_rh_medial_tmp.png'];
system(command);
command = ['freeview -f $CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/rh.very_inflated:annot=' ...
    rh_annot ':color=5,25,5:edgethickness=0 --viewport 3d -cam dolly 1.6 Azimuth 180 -ss ' ...
    fig_title '_rh_lateral_tmp.png'];
system(command);

%% trim the white bezel of the output images
cmd = CBIG_TrimFreeviewScreenshotCommand([fig_title '_lh_lateral_tmp.png'], ...
    [fig_save_path '/' fig_title '_lh_lateral.png']);
system(cmd);
cmd = CBIG_TrimFreeviewScreenshotCommand([fig_title '_lh_medial_tmp.png'], ...
    [fig_save_path '/' fig_title '_lh_medial.png']);
system(cmd);
cmd = CBIG_TrimFreeviewScreenshotCommand([fig_title '_rh_lateral_tmp.png'], ...
[fig_save_path '/' fig_title '_rh_lateral.png']);
system(cmd);
cmd = CBIG_TrimFreeviewScreenshotCommand([fig_title '_rh_medial_tmp.png'], ...
    [fig_save_path '/' fig_title '_rh_medial.png']);
system(cmd);

close all

delete([fig_title '_lh_lateral_tmp.png']);
delete([fig_title '_rh_lateral_tmp.png']);
delete([fig_title '_lh_medial_tmp.png']);
delete([fig_title '_rh_medial_tmp.png']);
delete('fsLR_164k_lh_annot.annot');
delete('fsLR_164k_rh_annot.annot');

end

