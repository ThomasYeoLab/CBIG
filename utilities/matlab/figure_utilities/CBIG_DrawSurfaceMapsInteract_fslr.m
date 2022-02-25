function CBIG_DrawSurfaceMapsInteract_fslr(lh_data, rh_data, mesh_name, surf_type, min_thresh, max_thresh, colors)

% CBIG_DrawSurfaceMapsInteract_fslr(lh_data, rh_data, mesh_name, surf_type, min_thresh, max_thresh, colors)
%
% This function is used to visualize a parcellation in fslr surface space but allow interactive
% cusor mode. By clicking a vertex in the figure, cusor will show the 
% vertex number and data of this vertex, e.g. if lh_data(9753)=0.5, then
% the cursor will show "9753:0.5"
%
% Input:
%      -lh_data, rh_data: 
%       surface labels on the left and right hemisphere respectively in
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
% Example:
% CBIG_DrawSurfaceMapsInteract_fslr(lh_labels,rh_labels,'fs_LR_32k','inflated', 0,400,rand(401,3));
%
% Written by Zhang Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% This function does not need vector check because the function itself
% contains checking statement.

warning('off', 'MATLAB:warn_r14_stucture_assignment');

global super_handles

if(~exist('mesh_name', 'var'))
   mesh_name = 'fs_LR_32k'; 
end

if(~exist('surf_type', 'var'))
   surf_type = 'inflated'; 
end

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


for hemis = {'lh' 'rh'}
    
    hemi = hemis{1};
    mesh= CBIG_read_fslr_surface(hemi, mesh_name, surf_type, 'medialwall.annot');
    non_cortex = find(mesh.MARS_label == 1);  
    
    if(strcmp(hemi, 'lh'))
        data = single(lh_data);
        super_handles.lh_mesh = mesh;
    elseif(strcmp(hemi, 'rh'))
        data = single(rh_data);
        super_handles.rh_mesh = mesh;
    end
    
    % convert to row vector
    if(size(data, 1) ~= 1)
       data = data';  
    end
    
    % check data dimensions
    if (max(size(mesh.vertices)) ~= max(size(data)))
        error('data dimension is different between label and data')
    end
    
    % threshold
    data(data < min_thresh) = min_thresh;
    data(data > max_thresh) = max_thresh;
    data(non_cortex(1)) = min_thresh;
    data(non_cortex(2)) = max_thresh;


    if(strcmp(hemi, 'lh'))
        super_handles.lh_data = data;
    elseif(strcmp(hemi, 'rh'))
        super_handles.rh_data = data;
    end
    
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

% interaction mode
h = datacursormode(gcf);
set(h,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on');
datacursormode on

end


function txt = myupdatefcn(empt, event_obj)

global super_handles

h = get(gcf, 'Children');

for i = 2:5
   if(~isempty(find(findall(h(i)) == get(event_obj, 'Target'), 1)))
       % right hemisphere
       vertex = MARS_findPoint(super_handles.rh_mesh.vertices, get(event_obj,'Position')');
       txt = [num2str(vertex) ': ' num2str(super_handles.rh_data(vertex))];
       return;
   end
end

for i = 6:9
   if(~isempty(find(findall(h(i)) == get(event_obj, 'Target'), 1)))
       % left hemisphere
       vertex = MARS_findPoint(super_handles.lh_mesh.vertices, get(event_obj,'Position')');
       txt = [num2str(vertex) ': ' num2str(super_handles.lh_data(vertex))];
       return;
   end
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
        ncd(x, y, 1: 3) = m(floor(idxf(x, y)) + 1, :);
    end
end
end