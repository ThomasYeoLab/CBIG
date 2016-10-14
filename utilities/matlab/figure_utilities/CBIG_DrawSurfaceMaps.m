function h = CBIG_DrawSurfaceMaps(lh_data, rh_data, mesh_name, surf_type, min_thresh, max_thresh, colors)

% h = CBIG_DrawSurfaceMaps(lh_data, rh_data, mesh_name, surf_type, min_thresh, max_thresh, colors)
%
% This function visualizes a given surface data in freesurfer space. 
% Threshold can be defined by min_thresh and max_thresh.
%
% Input:
%      -lh_data, rh_data: 
%       data of left/right hemisphere. Nx1 vector for each, 
%       N = # of vertices in mesh_name.
%
%      -mesh_name:
%       Freesurfer mesh structure. For example, 'fsaverage5'.
%
%      -surf_type:
%       Freesurfer surface template. Can be 'inflated', 'sphere', or
%       'white'.
%
%      -min_thresh, max_thresh:
%       min and max threshold for lh_data and rh_data. If they are not
%       given, then there is no threshold.
%
%      -colors:
%       color map for visualizetion. A Lx3 matrix, where L is the number of
%       different colors for lh_data and rh_data. Each row is the R, G, B
%       value. If colors is not given, visualization color will be defined
%       by default Matlab colormap.
%
% Output:
%       -h:
%        output figure number.
%
% Example:
% CBIG_DrawSurfaceMaps(lh_labels,rh_labels,'fsaverage5','inflated');
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(~exist('mesh_name', 'var'))
   mesh_name = 'fsaverage'; 
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
gpos(3) = 1200; gpos(4) = 600; set(h, 'Position', gpos);
for hemis = {'lh' 'rh'}
    
    hemi = hemis{1};
    mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, surf_type, 'cortex');
    non_cortex = find(mesh.MARS_label == 1);  
    
    if(strcmp(hemi, 'lh'))
        data = single(lh_data);
    elseif(strcmp(hemi, 'rh'))
        data = single(rh_data);
    end
    
    % convert to row vector
    if(size(data, 1) ~= 1)
       data = data';  
    end
    
    % resample
    if(size(mesh.vertices, 2) ~= length(data)) % need to resample!
        if(length(data) == 10242)
            from_mesh = CBIG_ReadNCAvgMesh(hemi, 'fsaverage5', 'sphere', 'cortex');
            target_mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, 'sphere', 'cortex');
            data = MARS_linearInterpolate(target_mesh.vertices, from_mesh, data);
        else
            error(['Not handling ' num2str(length(data)) ' vertices']);
        end
    end
    
    % threshold
    if(exist('min_thresh', 'var'))
       data(data < min_thresh) = min_thresh;
       data(data > max_thresh) = max_thresh;
       data(non_cortex(1)) = min_thresh;
       data(non_cortex(2)) = max_thresh;
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
