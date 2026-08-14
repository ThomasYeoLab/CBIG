function CBIG_EIAD_draw_surface_map(parameter, roi_list, fig_title)
% CBIG_EIAD_DRAW_SURFACE_MAP  Visualize a per-region value on the fsaverage5 surface.
%
% Renders a Desikan-parcellation value map (e.g. regional slopes or mediation
% effects) on inflated fsaverage5 surfaces, using the 'cool' colormap.
%
% Inputs:
%   parameter : M x 1 vector of values to display, one per shown region
%               (M = number of ones in roi_list).
%   roi_list  : 72 x 1 binary vector selecting which Desikan regions to display.
%               Set a region to 0 to hide it (e.g. it did not survive FDR).
%               Regions 1, 5, 37 and 41 (the medial wall) must always be 0.
%   fig_title : title string for the figure.
%
% Requires the CBIG repository (https://github.com/ThomasYeoLab/CBIG) on the
% MATLAB path (CBIG_ReadNCAvgMesh, TrisurfMeshData, MARS_*Interpolate, ...) and
% the file '1000subjects_clusters007_ref.mat' (which provides lh_labels, rh_labels).
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Map the displayed values onto the full 72-region vector
populated_vector = nan(72, 1);
counter = 1;
for i = 1:72
   if roi_list(i) == 1 
      populated_vector(i) = parameter(counter);
      counter = counter + 1;
   end
end
populated_vector([1, 5, 37, 41]) = [];   % drop the medial-wall regions

colorname = 'cool';

%% Load Desikan labels on fsaverage5
lh_mesh_fsavg = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'aparc.annot');
rh_mesh_fsavg = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'aparc.annot');

lh_label_desikan = lh_mesh_fsavg.MARS_label;
rh_label_desikan = rh_mesh_fsavg.MARS_label;


%% Load reference labels (provides lh_labels, rh_labels)
load('1000subjects_clusters007_ref.mat');

lh_label_desikan_data = lh_label_desikan;
rh_label_desikan_data = rh_label_desikan + 36;

for i = 1:72
   if i < 37
       if roi_list(i) == 0
           lh_label_desikan_data(lh_label_desikan_data == i) = 0;
       end
   else
       if roi_list(i) == 0
           rh_label_desikan_data(rh_label_desikan_data == i) = 0;
       end

   end
end

%% Build the per-region color table
HW = CBIG_EIAD_Num2Color(parameter,colorname);
HW = squeeze(HW(:,1,:));
HW = HW.*255;
index = 1:1:72;
index(~roi_list) = [];
HW_1 = ones(72,3);
HW_1(index,:) = HW;
HW_1 = [169 169 169; HW_1; 1 1 1];

%draw
visualize_surface(lh_label_desikan_data,rh_label_desikan_data, lh_labels, rh_labels, ...
    'fsaverage5', 'inflated', 0, 73, HW_1, populated_vector, fig_title)
end

function H = CBIG_EIAD_Num2Color(Wvector,colormap_name)
% Map a numeric vector to RGB colors using the (flipped) named colormap.
figure
colormap(colormap_name);
C = flip(colormap);
L = size(C,1);
[~, Windex] = sort(Wvector);
H_max = Wvector(Windex(end));
H_min = Wvector(Windex(1));
Ws = round(interp1(linspace(H_min,H_max,L),1:L,Wvector));
H = reshape(C(Ws,:),[size(Ws) 3]);
close
end

function visualize_surface(lh_data, rh_data, ...
    lh_labels, rh_labels, mesh_name, surf_type, min_thresh, max_thresh, colors, parameter, fig_title)

% adapted from CBIG_DrawSurfaceMapsWithBoundary

warning('off', 'MATLAB:warn_r14_stucture_assignment');

if(~exist('mesh_name', 'var'))
   mesh_name = 'fsaverage'; 
end

if(~exist('surf_type', 'var'))
   surf_type = 'inflated'; 
end

pos = [0.020, 0.510, 0.415, 0.470;...
    0.455, 0.510, 0.415, 0.470;...
    0.720, 0.760, 0.240, 0.230;...
    0.720, 0.510, 0.240, 0.230;...
    0.020, 0.020, 0.415, 0.470;...
    0.455, 0.020, 0.415, 0.470;...
    0.720, 0.260, 0.240, 0.230;...
    0.720, 0.010, 0.240, 0.230];

h = figure; gpos = get(h, 'Position');
gpos(1) = 100; gpos(2) = 100; gpos(3) = 800; gpos(4) = 600; set(h, 'Position', gpos);

if(exist('colors', 'var'))
    m = colors/max(colors(:));
    colormap(m);
else
    m = jet;
    colormap(m);
end

for hemis = {'lh' 'rh'}
    
    hemi = hemis{1};
    mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, surf_type, 'cortex');
    non_cortex = mesh.MARS_label == 1;  
    
    if(strcmp(hemi, 'lh'))
        data   = lh_data;
        labels = lh_labels;
    elseif(strcmp(hemi, 'rh'))
        data   = rh_data;
        labels = rh_labels;
    end

    % convert to row vector
    if(size(data, 1) ~= 1)
       data = data';  
    end
    
    if(size(labels, 1) ~= 1)
       labels = labels';  
    end
    
    % resample
    if(size(mesh.vertices, 2) ~= length(data)) % need to resample!
        if(length(data) == 10242)
            from_mesh = CBIG_ReadNCAvgMesh(hemi, 'fsaverage5', 'sphere', 'cortex');
            target_mesh = CBIG_ReadNCAvgMesh(hemi, mesh_name, 'sphere', 'cortex');
            data   = MARS_linearInterpolate(target_mesh.vertices, from_mesh, data);
            labels = MARS_NNInterpolate(target_mesh.vertices, from_mesh, labels);
        else
            error(['Not handling ' num2str(length(data)) ' vertices']);
        end
    end
    
    % threshold
    if(exist('min_thresh', 'var'))
       data(data < min_thresh) = min_thresh;
       data(data > max_thresh) = max_thresh;
       data(non_cortex) = min_thresh;
    end
    
    % compute boundary
    BoundaryVec = zeros(length(labels), 1);
    maxNeighbors = size(mesh.vertexNbors, 1);
    for i = 1:length(labels)
        label_vertex = int32(labels(i));
        
        for k = 1:maxNeighbors
            v_neighbor = mesh.vertexNbors(k, i);
            if(v_neighbor ~= 0 && int32(labels(v_neighbor)) ~= label_vertex)
                BoundaryVec(i) = 1;
            end
        end
    end
    data(BoundaryVec == 1) = max_thresh;
    
    % draw
    if(strcmp(hemi, 'lh'))
        subplot('Position', pos(1, :)); 
        s = TrisurfMeshData(mesh, data);
        shading('FLAT')
        view(-90, 0);
        axis off; 
        title(fig_title, 'Units', 'normalized', 'Position', [1, 1, 0], 'FontSize', 16)
        
        subplot('Position', pos(2, :));
        s = TrisurfMeshData(mesh, data);
        shading('FLAT')
        view(90, 0);
        axis off;
    else

        subplot('Position', pos(5, :));
        s = TrisurfMeshData(mesh, data);
        shading('FLAT')
        view(90, 0);
        axis off;

        subplot('Position', pos(6, :));
        s = TrisurfMeshData(mesh, data);
        shading('FLAT')
        view(-90, 0);
        axis off;
    end
end

ax1 = axes;
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
cmap = flipud(cool);
colormap(ax1, cmap);
set(gca, 'CLim', [min(parameter), max(parameter)])
lo = round(min(parameter), 3);
hi = round(max(parameter), 3);

colorbar(ax1, 'horiz', 'Position', [0.28 0.5 0.2 0.02], ...
    'XTick', [min(parameter), max(parameter)], ...
    'XTickLabel', {num2str(lo, '%g'), num2str(hi, '%g')});
end
