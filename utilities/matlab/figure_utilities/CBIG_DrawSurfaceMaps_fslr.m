function h = CBIG_DrawSurfaceMaps_fslr(lh_data, rh_data, mesh_name, surf_type, min_thresh,max_thresh,colors)

% h = CBIG_DrawSurfaceMaps_fslr(lh_data, rh_data, mesh_name, surf_type, min_thresh,max_thresh,colors)
%
% This function a given labels (parcellation labels, activation patterns...)
% in fslr to a Matlab figures of surfaces.
%
% Input:
%      - lh_data, rh_data:
%        surface labels on the left and right hemisphere respectively in
%        fslr. Each variable is a Nx1 or 1xN vector, where N corresponds to the
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

% This function does not need vector check because the function itself
% contains checking statement.

warning('off', 'MATLAB:warn_r14_stucture_assignment');

pos = [0.020, 0.510, 0.325, 0.470;...
    0.355, 0.510, 0.325, 0.470;...
    0.720, 0.760, 0.240, 0.230;...
    0.720, 0.510, 0.240, 0.230;...
    0.020, 0.020, 0.325, 0.470;...
    0.355, 0.020, 0.325, 0.470;...
    0.720, 0.260, 0.240, 0.230;...
    0.720, 0.010, 0.240, 0.230];

h = figure; gpos = get(h, 'Position');
gpos(1) = 0; gpos(2) = 0; gpos(3) = 1200; gpos(4) = 600; set(h, 'Position', gpos);

if(exist('colors', 'var'))
    m = colors/max(colors(:));
    colormap(m);
else
    m = jet;
    colormap(m);
end

%Add threshold if not specified
if(~exist('min_thresh', 'var'))
    min_thresh=min([min(lh_data) min(rh_data)]);
    max_thresh=max([max(lh_data) max(rh_data)]);
end


for hemis = {'lh','rh'}
    
    hemi = hemis{1};
    mesh= CBIG_read_fslr_surface(hemi, mesh_name, surf_type, 'medialwall.annot');
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
    
    %threshold
    data(data < min_thresh) = min_thresh;
    data(data > max_thresh) = max_thresh;
    data(non_cortex(1)) = min_thresh;
    data(non_cortex(2)) = max_thresh;
    
    
    % draw
    if(strcmp(hemi, 'lh'))
        subplot('Position', pos(1, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(-90, 0);
        axis off;
        
        subplot('Position', pos(2, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(90, 0);
        axis off;
        
        subplot('Position', pos(3, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(90, 90);
        axis off;
        
        subplot('Position', pos(8, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(90, -90);
        axis off;
        
    else
        
        subplot('Position', pos(5, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(90, 0);
        axis off;
        
        subplot('Position', pos(6, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(-90, 0);
        axis off;
        
        subplot('Position', pos(4, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(90, 90);
        axis off;
        
        subplot('Position', pos(7, :));
        s = TrisurfMeshData(mesh, data);
        shading interp;
        ncd = revert_shading_interp_behaviour(s, m);
        s.CData = ncd;
        view(90, -90);
        axis off;
    end
end

if(exist('min_thresh', 'var'))
    cbax = axes('Position', [0.29 0.5 0.1 0.02], 'visible', 'off');
    caxis(cbax, [min_thresh, max_thresh]);
    colorbar('peer', cbax, 'horiz', 'Position', [0.29 0.5 0.1 0.02]);
end

end


function ncd = revert_shading_interp_behaviour(s, m)
% shading interp behaviour is different across matlab versions
% we revert the shading interp behaviour to be like r2014a

s = get(s);
cdat = s.FaceVertexCData;
cl = get(gca, 'CLim');
sz = cl(2) - cl(1);
idxf = zeros(length(cdat), 1);
ncd = zeros(length(cdat), 1, 3);

for x = 1: length(cdat)
    for y = 1
        c = cdat(x, y);
        idxf(x, y) = ((c - cl(1)) / sz) * (size(m, 1) - 1);
        ncd(x, y, 1: 3) = m(round(idxf(x, y)) + 1, :);
    end
end
end

