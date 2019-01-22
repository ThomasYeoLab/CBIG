function h = CBIG_DrawSurfaceMaps_fslr(lh_data, rh_data, mesh_name, surf_type, min_thresh,max_thresh,colors)

% h = CBIG_DrawSurfaceMaps_fslr(lh_data, rh_data, mesh_name, surf_type, min_thresh,max_thresh,colors)
%
% This function a given labels (parcellation labels, activation patterns...)
% in fslr to a Matlab figures of surfaces.
%
% Input:
%      - lh_data, rh_data:
%        surface labels on the left and right hemisphere respectively in
%        fslr. Each variable is a Nx1 vector, where N corresponds to the
%        mesh defined in mesh_name.
%        N = 163,842 for mesh_name = 'fs_LR_164k'
%        N = 32,492 for mesh_name = 'fs_LR_32k'
%
%      - mesh_name:
%        fslr mesh that the label data will be projected onto.
%        fslr = 'fs_LR_32k' or 'fs_LR_164k'
%      - surf_type:
%        surface template. Options are 'inflated', 'very_inflated', 
%        'midthickness_orig', 'white_orig', 'pial_orig', 'sphere'.
%      - min_thresh, max_thresh (optional):
%        minimum and maximum threshold for the colorscale of the projected
%        values.
%      - colors (optional):
%        custom colorscale of the projected values
%
% Output:
%      - h:
%        handle of the figure
%
% Example:
% - CBIG_DrawSurfaceMaps_fslr(lh_proj_32K,rh_proj_32K, 'fs_LR_32k', 'inflated', 1e-05, 5e-5);
%   Draw label data onto fslr_32k inflated surface mesh with a colorscale
%   between 1e-05 and 5e-05.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


pos = [0.1 0.58 0.16 0.34; ...
    0.4 0.58 0.16 0.34; ...
    0.7 0.80 0.16 0.14; ...
    0.7 0.58 0.16 0.14; ...
    0.1 0.11 0.16 0.34; ...
    0.4 0.11 0.16 0.34; ...
    0.7 0.33 0.16 0.14; ...
    0.7 0.11 0.16 0.14];

h = figure; gpos = get(h, 'Position');
gpos(3) = 1200; gpos(4) = 600; set(h, 'Position', gpos);
for hemis = {'lh','rh'}
    
    hemi = hemis{1};
    mesh= CBIG_read_fslr_surface(hemi, mesh_name, surf_type);
  
    
    if(strcmp(hemi, 'lh'))
        data = single(lh_data);
    elseif(strcmp(hemi, 'rh'))
        data = single(rh_data);
    end
    
    non_cortex = find(data==0);  

    % convert to row vector
    if(size(data, 1) ~= 1)
       data = data';  
    end
    
    %threshold
    if(exist('min_thresh', 'var'))
       data(data < min_thresh) = min_thresh;
       data(data > max_thresh) = max_thresh;
       if numel(non_cortex) >= 2
           data(non_cortex(1)) = min_thresh;
           data(non_cortex(2)) = max_thresh;
       end
    end

    % draw
    if(strcmp(hemi, 'lh'))
        subplot('Position', pos(1, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(-90, 0); axis off; zoom(1.85);
        subplot('Position', pos(2, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(90, 0); axis off; zoom(1.85);
        subplot('Position', pos(3, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(90, 90); axis off; zoom(3.3);
        subplot('Position', pos(8, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(90, -90); axis off; zoom(3.3);
    else
        subplot('Position', pos(5, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(90, 0); axis off; zoom(1.85);
        subplot('Position', pos(6, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(-90, 0); axis off; zoom(1.85);
        subplot('Position', pos(4, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(90, 90); axis off; zoom(3.3);
        subplot('Position', pos(7, :)); TrisurfMeshData(mesh, data); shading interp; 
        view(90, -90); axis off; zoom(3.3);
    end
end

if(exist('colors', 'var'))
    colormap(colors/max(colors(:))); 
end

if(exist('min_thresh', 'var'))
    colorbar('horiz', 'Position', [0.28 0.5 0.1 0.02]);
end
