function CBIG_DrawfcMRISurfaceClusters(clusters_mat, mesh_name, surf_type, colors)

% CBIG_DrawfcMRISurfaceClusters(clusters_mat, mesh_name, surf_type, colors)
%
% This function is used to draw surface map for a given parcellation
%
% Input:
%      -clusters_mat:
%       a .mat file contains the parcellation information, include
%       lh_labels, rh_labels, colors
%
%      -mesh_name: 
%       mesh structure, e.g. fsaverage5
%
%      -surf_type: 
%       surface type of the template, can be 'inflated', 'sphere', 'white'
%
%      -colors: 
%       if not specified, the color map will be read from clusters_mat
%
% Example:
% CBIG_DrawfcMRISurfaceClusters(clusters, 'fsaverage5', 'inflated')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


clusters = load(clusters_mat);
if(~exist('colors', 'var')) && isfield(clusters, 'colors')
   colors = clusters.colors; 
end

if(~exist('mesh_name', 'var'))
   mesh_name = 'fsaverage5'; 
end

if(~exist('surf_type', 'var'))
   surf_type = 'inflated'; 
end

pos = [0.1 0.58 0.16 0.34; ...
    0.4 0.58 0.16 0.34; ...
    0.7 0.80 0.16 0.14; ...
    0.7 0.58 0.16 0.14; ...
    0.1 0.11 0.16 0.34; ...
    0.4 0.11 0.16 0.34; ...
    0.7 0.33 0.16 0.14; ...
    0.7 0.11 0.16 0.14];

h = figure; gpos = get(h, 'Position');
gpos(1) = 0; gpos(2) = 0; gpos(3) = 1200; gpos(4) = 600; set(h, 'Position', gpos);
for hemis = {'lh' 'rh'}
    
    hemi = hemis{1};
    mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, surf_type, 'cortex');
    
    if(strcmp(hemi, 'lh'))
        labels = clusters.lh_labels;
    elseif(strcmp(hemi, 'rh'))
        labels = clusters.rh_labels;
    end
    
    if(strcmp(hemi, 'lh'))
        subplot('Position', pos(1, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(-90, 0); axis off; zoom(1.85);
        subplot('Position', pos(2, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(90, 0); axis off; zoom(1.85);
        subplot('Position', pos(3, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(90, 90); axis off; zoom(3.3);
        subplot('Position', pos(8, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(90, -90); axis off; zoom(3.3);
    else
        subplot('Position', pos(5, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(90, 0); axis off; zoom(1.85);
        subplot('Position', pos(6, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(-90, 0); axis off; zoom(1.85);
        subplot('Position', pos(4, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(90, 90); axis off; zoom(3.3);
        subplot('Position', pos(7, :)); TrisurfMeshData(mesh, labels); shading flat; 
        view(90, -90); axis off; zoom(3.3);
    end
end

if(exist('colors', 'var'))
    colormap(colors/max(colors(:))); 
end