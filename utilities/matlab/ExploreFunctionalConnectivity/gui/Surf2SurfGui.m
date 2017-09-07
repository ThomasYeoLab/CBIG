function varargout = Surf2SurfGui(varargin)
% SURF2SURFGUI M-file for Surf2SurfGui.fig
%      SURF2SURFGUI, by itself, creates a new SURF2SURFGUI or raises the
%      existing
%      singleton*.
%
%      H = SURF2SURFGUI returns the handle to a new SURF2SURFGUI or the handle to
%      the existing singleton*.
%
%      SURF2SURFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SURF2SURFGUI.M with the given input arguments.
%
%      SURF2SURFGUI('Property','Value',...) creates a new SURF2SURFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Surf2SurfGui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Surf2SurfGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Edit the above text to modify the response to help Surf2SurfGui

% Last Modified by GUIDE v2.5 01-Nov-2010 18:09:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Surf2SurfGui_OpeningFcn, ...
                   'gui_OutputFcn',  @Surf2SurfGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Surf2SurfGui is made visible.
function Surf2SurfGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Surf2SurfGui (see VARARGIN)

% Choose default command line output for Surf2SurfGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

% add optional corr_mat input
[a, hostname] = system('hostname -d');
corr_mat_file = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'surf2surf_gui_data', 'PrecomputedSurfaceCorrelation', 'full_surface_corr.1000sub.mat');
p.addParamValue('corr_mat_file', corr_mat_file);

% add option for opengl or zbuffer renderer 
p.addParamValue('renderer', 'zbuffer', @(x)strcmp(x,'opengl') || ...
              strcmp(x,'zbuffer'));
          
% add optional number of networks. 0 ==> don't load networks
p.addParamValue('num_networks', 0);

% add dataset_num used to estimate the networks.
p.addParamValue('dataset_num', 3);

% optional input for silhouette
p.addParamValue('sil_mat_file', []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provide help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(sum(strcmp(varargin, 'help')) > 0)
   disp(['Surf2SurfGui <corr_mat_file corr_mat_file> <renderer renderer_string> <num_networks num_networks> <dataset_num dataset_num>']);
   disp(' ');
    disp('Example Calls: ');
    disp(['    Surf2SurfGui renderer opengl']);
    disp(['    Surf2SurfGui num_networks 7 corr_mat_file ' corr_mat_file]);
    disp(['    Surf2SurfGui renderer zbuffer dataset_num 3 num_networks 7 corr_mat_file ' corr_mat_file]);
    disp(' ');
    disp(['default corr_mat_file = ' corr_mat_file]);
    disp(['default renderer = zbuffer (alternative is opengl)']);
    disp('default dataset_num = 3, which uses parcellation results from 1000 subjects (data_num = 1 uses 500 subjects, data_num = 2 uses second 500 subjects)');
    disp('default number_networks = 0, which means use low resolution surface version without showing networks. Available choices = 7, 10, 12, 17');
    handles.help = 1;
    guidata(hObject, handles);
    return; 
else
    handles.help = 0;
    guidata(hObject, handles);
end

p.parse(varargin{:});
corr_mat_file = p.Results.corr_mat_file;
num_networks = p.Results.num_networks;
dataset_num = p.Results.dataset_num;
handles.renderer = p.Results.renderer;

if(ischar(num_networks))
    num_networks = str2num(num_networks);
end

if(ischar(dataset_num))
    dataset_num = str2num(dataset_num);
end

handles.corr_mat_file = corr_mat_file;
handles.num_networks = num_networks;
handles.dataset_num = dataset_num;
handles.sil_mat_file = p.Results.sil_mat_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin data reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read mesh
disp('Reading meshes ...');
if(num_networks > 0)
    handles.lh_avg_mesh5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'sphere', 'cortex');
    handles.rh_avg_mesh5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'sphere', 'cortex'); 
    handles.lh_avg_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
    handles.rh_avg_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex'); 
    handles.lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'cortex');
    handles.rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'cortex'); 
    handles.lh_avg_mesh_inflated5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
    handles.rh_avg_mesh_inflated5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex'); 
else
    handles.lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
    handles.rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex'); 
    handles.lh_avg_mesh_inflated5 = handles.lh_avg_mesh;
    handles.rh_avg_mesh_inflated5 = handles.rh_avg_mesh;
end
handles.curv_gray_scale = [0.65 0.65 0.65; 0.4 0.4 0.4; 0 0 0];
low_threshold = str2double(get(handles.low_threshold,'String'));
high_threshold = str2double(get(handles.high_threshold,'String'));

% binarizing curvature
lh_curv = handles.lh_avg_mesh.data(2, :);
lh_curv(lh_curv >= 0) = 1.5;
lh_curv(lh_curv < 0) = 1.25;
lh_curv(handles.lh_avg_mesh.MARS_label == 1) = 2;
lh_curv(1) = 0;
handles.lh_curv = lh_curv;
    
rh_curv = handles.rh_avg_mesh.data(2, :);
rh_curv(rh_curv >= 0) = 1.5;
rh_curv(rh_curv < 0) = 1.25;
rh_curv(handles.rh_avg_mesh.MARS_label == 1) = 2;
rh_curv(1) = 0;
handles.rh_curv = rh_curv;

% read parcellations
if(num_networks > 0)
    disp('Reading parcellations');
    num_networks_str = num2str(num_networks, '%03d');
    num_networks_str_prev = num2str(num_networks-1, '%03d');

    [a, hostname] = system('hostname -d');
    data_path = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'surf2surf_gui_data', 'cluster_results', ['cluster' num_networks_str]);
    colortable_path = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'figure_utilities');

    if(dataset_num == 2)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.500sub.2.ref500sub.1.mat']);
    elseif(dataset_num == 1)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '_ref' num_networks_str_prev '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.500sub.1.mat']);
    elseif(dataset_num == 3)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.1000sub.ref500sub.1.mat']);
    else
        error(['Does not handle data set num = ' num2str(dataset_num)]);
    end
    load(network_file);
    handles.lh_labels = lh_labels';
    handles.rh_labels = rh_labels';
    
    % upsample labels to high resolution
    handles.lh_labels7 = MARS_NNInterpolate(handles.lh_avg_mesh7.vertices, handles.lh_avg_mesh5, handles.lh_labels);
    handles.rh_labels7 = MARS_NNInterpolate(handles.rh_avg_mesh7.vertices, handles.rh_avg_mesh5, handles.rh_labels);

    load(fullfile(colortable_path, 'freesurfer_color.mat'));
    handles.freesurfer_color = freesurfer_color(1:num_networks+1, :);

    % silhouettte
    if(dataset_num == 2)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.500sub.2.mat']);
    elseif(dataset_num == 1)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.500sub.1.mat']);
    elseif(dataset_num == 3)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.1000sub.mat']);
    else
        error(['Does not handle data set num = ' num2str(dataset_num)]);
    end
    load(network_file);
    
    if(~isempty(handles.sil_mat_file))
       disp(['Reading silhouette file: ' handles.sil_mat_file]);
       load(handles.sil_mat_file); 
    end
    lh_s(lh_s == 1) = min(lh_s);
    rh_s(rh_s == 1) = min(rh_s);
    handles.lh_sil = lh_s';
    handles.rh_sil = rh_s';
    
    % upsample sil to high resolution and precompute interpolation weights
    [handles.lh_sil7, tmp, NF] = MARS_linearInterpolate(handles.lh_avg_mesh7.vertices, handles.lh_avg_mesh5, handles.lh_sil);
    [handles.lh_weights, handles.lh_vno] = MARS_computeWeightsForInterpolation(handles.lh_avg_mesh5, handles.lh_avg_mesh7.vertices, NF);
    [handles.rh_sil7, tmp, NF] = MARS_linearInterpolate(handles.rh_avg_mesh7.vertices, handles.rh_avg_mesh5, handles.rh_sil);
    [handles.rh_weights, handles.rh_vno] = MARS_computeWeightsForInterpolation(handles.rh_avg_mesh5, handles.rh_avg_mesh7.vertices, NF);
    
    % precompute boundary vec
    maxNeighbors = size(handles.lh_avg_mesh7.vertexNbors, 1);
    handles.lh_boundary_vec = false(1, length(handles.lh_labels7));
    for i = 1:length(handles.lh_labels7)
        label_vertex = int32(handles.lh_labels7(i));

        for k = 1:maxNeighbors
            v_neighbor = handles.lh_avg_mesh7.vertexNbors(k, i);
            if(v_neighbor ~= 0 && int32(handles.lh_labels7(v_neighbor)) ~= label_vertex)
                handles.lh_boundary_vec(i) = 1;
            end
        end
    end

    maxNeighbors = size(handles.rh_avg_mesh7.vertexNbors, 1);
    handles.rh_boundary_vec = false(1, length(handles.rh_labels7));
    for i = 1:length(handles.rh_labels7)
        label_vertex = int32(handles.rh_labels7(i));

        for k = 1:maxNeighbors
            v_neighbor = handles.rh_avg_mesh7.vertexNbors(k, i);
            if(v_neighbor ~= 0 && int32(handles.rh_labels7(v_neighbor)) ~= label_vertex)
                handles.rh_boundary_vec(i) = 1;
            end
        end
    end
end
% create new colormap;
handles.corr_colormap = CreateCorrColormap(handles.curv_gray_scale, low_threshold, high_threshold);
    
% reading correlation matrix
disp(['Loading correlation mat: ' corr_mat_file]);
corr_mat = load(corr_mat_file, 'corr_mat');
handles.corr_mat = corr_mat.corr_mat;
clear corr_mat;

% reading coordinate systems
disp('Loading Coordinate Transforms');
coord_path = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Buckner2011');

load(fullfile(coord_path, 'coord_surf2vol', 'lh.1000sub.FS.vox.mat')); handles.lh_FS_vox = vox;
load(fullfile(coord_path, 'coord_surf2vol', 'rh.1000sub.FS.vox.mat')); handles.rh_FS_vox = vox;

load(fullfile(coord_path, 'coord_surf2vol', 'lh.1000sub.FS.ras.mat')); handles.lh_FS_ras = ras;
load(fullfile(coord_path, 'coord_surf2vol', 'rh.1000sub.FS.ras.mat')); handles.rh_FS_ras = ras;

load(fullfile(coord_path, 'coord_surf2vol', 'lh.1000sub.FSL_MNI.vox.mat')); handles.lh_FSL_MNI152_vox = vox;
load(fullfile(coord_path, 'coord_surf2vol', 'rh.1000sub.FSL_MNI.vox.mat')); handles.rh_FSL_MNI152_vox = vox;

load(fullfile(coord_path, 'coord_surf2vol', 'lh.1000sub.FSL_MNI.ras.mat')); handles.lh_FSL_MNI152_ras = ras;
load(fullfile(coord_path, 'coord_surf2vol', 'rh.1000sub.FSL_MNI.ras.mat')); handles.rh_FSL_MNI152_ras = ras;

handles.FS_dist = MRIread(fullfile(coord_path, 'coord_vol2surf', '1000sub.FS.2mm.dist_map.500.nii.gz'));
handles.FS_vertex = MRIread(fullfile(coord_path, 'coord_vol2surf', '1000sub.FS.2mm.full_vertex_map.500.fsaverage5.nii.gz'));
handles.FSL_MNI152_dist = MRIread(fullfile(coord_path, 'coord_vol2surf', '1000sub.FSL_MNI152.1mm.dist_map.500.nii.gz'));
handles.FSL_MNI152_vertex = MRIread(fullfile(coord_path, 'coord_vol2surf', '1000sub.FSL_MNI152.1mm.full_vertex_map.500.fsaverage5.nii.gz'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start drawing!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(num_networks > 0)
    max_sil = max(max(handles.lh_sil), max(handles.rh_sil));
    min_sil = min(min(handles.lh_sil), min(handles.rh_sil));
    jet_colormap = jet(128);
    num_high_colors = round(size(jet_colormap, 1) / (max_sil - min_sil) * (1 - max_sil) * 0.5); % 0.5
    num_low_colors  = round(size(jet_colormap, 1) / (max_sil - min_sil) * (min_sil + 1)); 
    padded_jet_colormap = [repmat(jet_colormap(1, :), num_low_colors, 1) ; jet_colormap; repmat(jet_colormap(end, :), num_high_colors, 1)];
    
    num_colors_per_label = round(size(jet_colormap, 1) / (max_sil - min_sil));
    num_remain_colors = 2*num_colors_per_label - size(padded_jet_colormap, 1);
    padded_colormap = [];
    for i = 2:num_networks
        padded_colormap = [padded_colormap; repmat(handles.freesurfer_color(i, :), num_colors_per_label, 1)];
    end
    num_last_network_colors = round(num_colors_per_label/2);
    num_first_network_colors = num_colors_per_label - num_last_network_colors;
    padded_colormap = [repmat(handles.freesurfer_color(1, :), num_first_network_colors, 1); padded_colormap; repmat(handles.freesurfer_color(num_networks+1, :), num_last_network_colors, 1)];
    
    handles.label_sil_colormap = [padded_jet_colormap; repmat(handles.freesurfer_color(1, :), num_remain_colors, 1); padded_colormap];
    
    axes(handles.LeftHemisphere)
    tmp = handles.lh_labels+1;
    tmp(1) = -1;
    TrisurfMeshData(handles.lh_avg_mesh_inflated5, tmp);
    shading flat; colormap(handles.label_sil_colormap);
    axis off;
    view(-90, 0);
    
    axes(handles.RightHemisphere)
    tmp = handles.rh_labels+1;
    tmp(1) = -1;
    TrisurfMeshData(handles.rh_avg_mesh_inflated5, tmp);
    shading flat; colormap(handles.label_sil_colormap);
    axis off;
    view(90, 0);
    
    axes(handles.LeftHemisphereSil)
    tmp = handles.lh_sil;
    tmp(1) = -1;
    tmp(2) = num_networks+1;
    TrisurfMeshData(handles.lh_avg_mesh_inflated5, tmp);
    shading interp; colormap(handles.label_sil_colormap);
    axis off;
    view(-90, 0);
    
    axes(handles.RightHemisphereSil)
    tmp = handles.rh_sil;
    tmp(1) = -1;
    tmp(2) = num_networks+1;
    TrisurfMeshData(handles.rh_avg_mesh_inflated5, tmp);
    shading interp; colormap(handles.label_sil_colormap);
    axis off;
    view(90, 0);
else
    axes(handles.LeftHemisphere)
    TrisurfMeshData(handles.lh_avg_mesh_inflated5, handles.lh_curv);
    colormap(handles.corr_colormap);
    shading interp;
    axis off;
    view(-90, 0);
    
    axes(handles.RightHemisphere)
    TrisurfMeshData(handles.rh_avg_mesh_inflated5, handles.rh_curv);
    colormap(handles.corr_colormap);
    shading interp;
    axis off;
    view(90, 0);
    
    axes(handles.LeftHemisphereSil); axis off;
    axes(handles.RightHemisphereSil); axis off;
end


% Set up data cursor structure
dcm_obj = datacursormode(handles.figure1);
handles.dcm_obj = dcm_obj;
set(dcm_obj,'UpdateFcn',@myupdatefcn)

% start with rotate mode.
rotate3d on;

% initialize vertex
handles.vertex = GetVertexFromVoxel(str2num(get(handles.FS_vox, 'String')), handles.FS_vertex);
if(handles.vertex < 0)
   handles.vertex = 10242 + abs(vertex); 
end
UpdateFigures(handles, handles.vertex)

global super_handles; super_handles = handles;
guidata(hObject, handles);


function txt = myupdatefcn(empt, event_obj)
% Customizes text of data tips

global super_handles;

if(~isempty(find(findall(super_handles.LeftHemisphere) == get(event_obj, 'Target'), 1)) || ~isempty(find(findall(super_handles.LeftHemisphereSil) == get(event_obj, 'Target'), 1)))
    hemi = 'lh';
    vertex = MARS_findPoint(super_handles.lh_avg_mesh_inflated5.vertices, get(event_obj,'Position')');
    tmp = super_handles.lh_FS_ras(:, vertex);
    if(super_handles.num_networks > 0)
        txt = {['FS_X: ',num2str(tmp(1))],...
            ['FS_Y: ', num2str(tmp(2))], ...
            ['FS_Z: ', num2str(tmp(3))], ...
            ['network: ', num2str(super_handles.lh_labels(vertex))], ...
            ['sil: ', num2str(super_handles.lh_sil(vertex))]};
    else
        txt = {['FS_X: ',num2str(tmp(1))],...
            ['FS_Y: ', num2str(tmp(2))], ...
            ['FS_Z: ', num2str(tmp(3))]};

    end
elseif(~isempty(find(findall(super_handles.RightHemisphere) == get(event_obj, 'Target'), 1)) || ~isempty(find(findall(super_handles.RightHemisphereSil) == get(event_obj, 'Target'), 1)))
    hemi = 'rh';
    vertex = MARS_findPoint(super_handles.rh_avg_mesh_inflated5.vertices, get(event_obj,'Position')');
    tmp = super_handles.rh_FS_ras(:, vertex);
    if(super_handles.num_networks > 0)
        txt = {['FS_X: ',num2str(tmp(1))],...
            ['FS_Y: ', num2str(tmp(2))], ...
            ['FS_Z: ', num2str(tmp(3))], ...
            ['network: ', num2str(super_handles.rh_labels(vertex))], ...
            ['sil: ', num2str(super_handles.rh_sil(vertex))]};
    else
        txt = {['FS_X: ',num2str(tmp(1))],...
            ['FS_Y: ', num2str(tmp(2))], ...
            ['FS_Z: ', num2str(tmp(3))]};
    end
else
    warning('Please select point from either figures on the left or right');
end


% UIWAIT makes Surf2SurfGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function corr_colormap = CreateCorrColormap(curv_gray_scale, low_threshold, high_threshold)

jet_colormap = jet(64);
num_high_colors = round(size(jet_colormap, 1) / (high_threshold - low_threshold) * (1 - high_threshold));
num_low_colors  = round(size(jet_colormap, 1) / (high_threshold - low_threshold) * low_threshold);
padded_jet_colormap = [repmat(jet_colormap(1, :), num_low_colors, 1) ; jet_colormap; repmat(jet_colormap(end, :), num_high_colors, 1)];
num_colors = size(padded_jet_colormap, 1);

num_low_colors = round(num_colors/3);
num_mid_colors = round(num_colors/3);
num_high_colors = num_colors - num_low_colors - num_mid_colors;
corr_colormap = [padded_jet_colormap; repmat(curv_gray_scale(1, :), num_low_colors, 1); repmat(curv_gray_scale(2, :), num_mid_colors, 1); repmat(curv_gray_scale(3, :), num_high_colors, 1)];



% --- Outputs from this function are returned to the command line.
function varargout = Surf2SurfGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

if(handles.help)
    delete(handles.figure1);
end

function FS_vox_Callback(hObject, eventdata, handles)
% hObject    handle to FS_vox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FS_vox as text
%        str2double(get(hObject,'String')) returns contents of FS_vox as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

vox = round(str2num(get(hObject, 'String')));
vox = vox';
handles.vertex = GetVertexFromVoxel(vox, handles.FS_vertex);
if(handles.vertex < 0)
    set(handles.rh_vertex, 'String', num2str(abs(handles.vertex)));
    set(handles.lh_vertex, 'String', 'rh vertex');
    handles.vertex = 10242 + abs(handles.vertex);
else
    set(handles.lh_vertex, 'String', num2str(vertex));
    set(handles.rh_vertex, 'String', 'lh vertex');
end

ras = CBIG_ConvertVox2Ras(vox, handles.FS_vertex.vox2ras);
set(handles.FS_vox, 'String', [num2str(vox(1)) ', ' num2str(vox(2)) ', ' num2str(vox(3))]);
set(handles.FS_ras, 'String', [num2str(ras(1)) ', ' num2str(ras(2)) ', ' num2str(ras(3))]);

set(handles.FSL_MNI152_vox, 'String', 'FS mode');
set(handles.FSL_MNI152_ras, 'String', 'FS mode');

set(handles.FSL_MNI152_DistToCortex, 'String', ['Dist to cortex: FS mode']);
set(handles.FS_DistToCortex, 'String', ['Dist to cortex: ' num2str(handles.FS_dist.vol(vox(2)+1, vox(1)+1, vox(3)+1)) ' mm']);

UpdateFigures(handles, handles.vertex)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FS_vox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FS_vox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function FS_ras_Callback(hObject, eventdata, handles)
% hObject    handle to FS_ras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FS_ras as text
%        str2double(get(hObject,'String')) returns contents of FS_ras as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

% this is necessary because we are in 2mm space, so we want to round off to
% the nearest even number
ras = round(str2num(get(hObject, 'String'))/2)*2;
ras = ras';
handles.vertex = GetVertexFromRAS(ras, handles.FS_vertex);
if(handles.vertex < 0)
    set(handles.rh_vertex, 'String', num2str(abs(handles.vertex)));
    set(handles.lh_vertex, 'String', 'rh vertex');
    handles.vertex = 10242 + abs(handles.vertex);
else
    set(handles.lh_vertex, 'String', num2str(handles.vertex));
    set(handles.rh_vertex, 'String', 'lh vertex');
end

vox = CBIG_ConvertRas2Vox(ras, handles.FS_vertex.vox2ras);
set(handles.FS_vox, 'String', [num2str(vox(1)) ', ' num2str(vox(2)) ', ' num2str(vox(3))]);
set(handles.FS_ras, 'String', [num2str(ras(1)) ', ' num2str(ras(2)) ', ' num2str(ras(3))]);

set(handles.FSL_MNI152_vox, 'String', 'FS mode');
set(handles.FSL_MNI152_ras, 'String', 'FS mode');

set(handles.FSL_MNI152_DistToCortex, 'String', ['Dist to cortex: FS mode']);
set(handles.FS_DistToCortex, 'String', ['Dist to cortex: ' num2str(handles.FS_dist.vol(vox(2)+1, vox(1)+1, vox(3)+1)) ' mm']);

UpdateFigures(handles, handles.vertex)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FS_ras_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FS_ras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FSL_MNI152_vox_Callback(hObject, eventdata, handles)
% hObject    handle to FSL_MNI152_vox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FSL_MNI152_vox as text
%        str2double(get(hObject,'String')) returns contents of FSL_MNI152_vox as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

vox = round(str2num(get(hObject, 'String')));
vox = vox';
handles.vertex = GetVertexFromVoxel(vox, handles.FSL_MNI152_vertex);
if(handles.vertex < 0)
    set(handles.rh_vertex, 'String', num2str(abs(handles.vertex)));
    set(handles.lh_vertex, 'String', 'rh vertex');
    handles.vertex = 10242 + abs(handles.vertex);
else
    set(handles.lh_vertex, 'String', num2str(vertex));
    set(handles.rh_vertex, 'String', 'lh vertex');
end

ras = CBIG_ConvertVox2Ras(vox, handles.FSL_MNI152_vertex.vox2ras);

set(handles.FSL_MNI152_vox, 'String', [num2str(vox(1)) ', ' num2str(vox(2)) ', ' num2str(vox(3))]);
set(handles.FSL_MNI152_ras, 'String', [num2str(ras(1)) ', ' num2str(ras(2)) ', ' num2str(ras(3))]);

set(handles.FS_vox, 'String', 'FSL MNI152 mode');
set(handles.FS_ras, 'String', 'FSL MNI152 mode');

set(handles.FS_DistToCortex, 'String', ['Dist to cortex: FSL MNI152 mode']);
set(handles.FSL_MNI152_DistToCortex, 'String', ['Dist to cortex: ' num2str(handles.FSL_MNI152_dist.vol(vox(2)+1, vox(1)+1, vox(3)+1)) ' mm']);

UpdateFigures(handles, handles.vertex)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FSL_MNI152_vox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FSL_MNI152_vox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FSL_MNI152_ras_Callback(hObject, eventdata, handles)
% hObject    handle to FSL_MNI152_ras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FSL_MNI152_ras as text
%        str2double(get(hObject,'String')) returns contents of FSL_MNI152_ras as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

% this is necessary because we are in 2mm space, so we want to round off to
% the nearest even number
ras = round(str2num(get(hObject, 'String')));
ras = ras';
handles.vertex = GetVertexFromRAS(ras, handles.FSL_MNI152_vertex);
vertex = handles.vertex;
if(handles.vertex < 0)
    set(handles.rh_vertex, 'String', num2str(abs(handles.vertex)));
    set(handles.lh_vertex, 'String', 'rh vertex');
    handles.vertex = 10242 + abs(handles.vertex);
else
    set(handles.lh_vertex, 'String', num2str(vertex));
    set(handles.rh_vertex, 'String', 'lh vertex');
end

vox = CBIG_ConvertRas2Vox(ras, handles.FSL_MNI152_vertex.vox2ras);
set(handles.FSL_MNI152_vox, 'String', [num2str(vox(1)) ', ' num2str(vox(2)) ', ' num2str(vox(3))]);
set(handles.FSL_MNI152_ras, 'String', [num2str(ras(1)) ', ' num2str(ras(2)) ', ' num2str(ras(3))]);

set(handles.FS_vox, 'String', 'FSL_MNI152 mode');
set(handles.FS_ras, 'String', 'FSL_MNI152 mode');

set(handles.FS_DistToCortex, 'String', ['Dist to cortex: FSL_MNI152 mode']);
set(handles.FSL_MNI152_DistToCortex, 'String', ['Dist to cortex: ' num2str(handles.FSL_MNI152_dist.vol(vox(2)+1, vox(1)+1, vox(3)+1)) ' mm']);

UpdateFigures(handles, handles.vertex)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FSL_MNI152_ras_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FSL_MNI152_ras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function low_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to low_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_threshold as text
%        str2double(get(hObject,'String')) returns contents of low_threshold as a double

low_threshold = str2double(get(handles.low_threshold,'String'));
high_threshold = str2double(get(handles.high_threshold,'String'));

if(high_threshold < low_threshold)
   disp('high threshold should be greater than low threshold'); 
   set(handles.low_threshold,'String',get(handles.high_threshold,'String'));
end

handles.corr_colormap = CreateCorrColormap(handles.curv_gray_scale, low_threshold, high_threshold);
guidata(hObject, handles);
if(~isempty(handles.vertex))
    UpdateFigures(handles, handles.vertex);
end


% --- Executes during object creation, after setting all properties.
function low_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function high_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to high_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high_threshold as text
%        str2double(get(hObject,'String')) returns contents of high_threshold as a double

low_threshold = str2double(get(handles.low_threshold,'String'));
high_threshold = str2double(get(handles.high_threshold,'String'));

if(high_threshold < low_threshold)
   disp('high threshold should be greater than low threshold'); 
   set(handles.high_threshold,'String',get(handles.low_threshold,'String'));
end

handles.corr_colormap = CreateCorrColormap(handles.curv_gray_scale, low_threshold, high_threshold);
guidata(hObject, handles);
if(~isempty(handles.vertex))
    UpdateFigures(handles, handles.vertex);
end


% --- Executes during object creation, after setting all properties.
function high_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to high_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Rotate.
function Rotate_Callback(hObject, eventdata, handles)
% hObject    handle to Rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Rotate
set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
rotate3d on;

% --- Executes on button press in SelectVertex.
function SelectVertex_Callback(hObject, eventdata, handles)
% hObject    handle to SelectVertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SelectVertex
rotate3d off;
while(get(hObject, 'Value'))
  set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on');
  
  disp('Select vertex and press enter');  
  pause                            % Wait while the user does this.
  cursor_info = getCursorInfo(handles.dcm_obj);
  if(~isempty(cursor_info))
      handles = guidata(gcbo); % this is necessary because stuff like the colormap might have changed between clicking vertices 
      
      if(~isempty(find(findall(handles.LeftHemisphere) == cursor_info.Target, 1)) || ~isempty(find(findall(handles.LeftHemisphereSil) == cursor_info.Target, 1)))
          hemi = 'lh';
          vertex = MARS_findPoint(handles.lh_avg_mesh_inflated5.vertices, cursor_info.Position');
          set(handles.lh_vertex, 'String', num2str(vertex));
          set(handles.rh_vertex, 'String', 'lh vertex');
          
          % Update freesurfer vox, ras
          set(handles.FS_vox, 'String', [num2str(handles.lh_FS_vox(1, vertex), '%.1f') ', ' num2str(handles.lh_FS_vox(2, vertex), '%.1f') ', ' num2str(handles.lh_FS_vox(3, vertex), '%.1f')]);
          set(handles.FS_ras, 'String', [num2str(handles.lh_FS_ras(1, vertex), '%.1f') ', ' num2str(handles.lh_FS_ras(2, vertex), '%.1f') ', ' num2str(handles.lh_FS_ras(3, vertex), '%.1f')]);  
          
          % Update FSL_MNI152 vox, ras
          set(handles.FSL_MNI152_vox, 'String', [num2str(handles.lh_FSL_MNI152_vox(1, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_vox(2, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_vox(3, vertex), '%.1f')]);
          set(handles.FSL_MNI152_ras, 'String', [num2str(handles.lh_FSL_MNI152_ras(1, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_ras(2, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_ras(3, vertex), '%.1f')]);  
      elseif(~isempty(find(findall(handles.RightHemisphere) == cursor_info.Target, 1)) || ~isempty(find(findall(handles.RightHemisphereSil) == cursor_info.Target, 1)))
          hemi = 'rh';  
          vertex = MARS_findPoint(handles.rh_avg_mesh_inflated5.vertices, cursor_info.Position');
          set(handles.rh_vertex, 'String', num2str(vertex));
          set(handles.lh_vertex, 'String', 'rh vertex');
          
          % Update freesurfer vox, ras
          set(handles.FS_vox, 'String', [num2str(handles.rh_FS_vox(1, vertex), '%.1f') ', ' num2str(handles.rh_FS_vox(2, vertex), '%.1f') ', ' num2str(handles.rh_FS_vox(3, vertex), '%.1f')]);
          set(handles.FS_ras, 'String', [num2str(handles.rh_FS_ras(1, vertex), '%.1f') ', ' num2str(handles.rh_FS_ras(2, vertex), '%.1f') ', ' num2str(handles.rh_FS_ras(3, vertex), '%.1f')]);  
          
          % Update FSL_MNI152 vox, ras
          set(handles.FSL_MNI152_vox, 'String', [num2str(handles.rh_FSL_MNI152_vox(1, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_vox(2, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_vox(3, vertex), '%.1f')]);
          set(handles.FSL_MNI152_ras, 'String', [num2str(handles.rh_FSL_MNI152_ras(1, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_ras(2, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_ras(3, vertex), '%.1f')]);  
          vertex = vertex + 10242;
      else
         warning('Please select point from either figures on the left');
         continue
      end
      handles.vertex = vertex;
      guidata(hObject, handles);
      UpdateFigures(handles, vertex);
      
      % Update distance. 
      set(handles.FS_DistToCortex, 'String', 'Dist to cortex : surf mode');
      set(handles.FSL_MNI152_DistToCortex, 'String', 'Dist to cortex : surf mode');
      
      
  end
end

function UpdateFigures(handles, vertex)

h = figure(2); set(gcf, 'renderer', handles.renderer);
pos = get(h, 'Position');
pos(3) = 1000;
pos(4) = 700;
set(h, 'Position', pos);

lh_corr = handles.corr_mat(vertex, 1:10242);
rh_corr = handles.corr_mat(vertex, 10243:end);

low_threshold = str2double(get(handles.low_threshold,'String'));
high_threshold = str2double(get(handles.high_threshold,'String'));

lh_corr(handles.lh_avg_mesh_inflated5.MARS_label == 1) = 0; lh_corr(lh_corr < low_threshold) = 0; lh_corr(lh_corr > high_threshold) = high_threshold;
rh_corr(handles.rh_avg_mesh_inflated5.MARS_label == 1) = 0; rh_corr(rh_corr < low_threshold) = 0; rh_corr(rh_corr > high_threshold) = high_threshold;

fprintf('Please wait ...')
if(handles.num_networks > 0)
    lh_corr = sum(lh_corr(handles.lh_vno) .* handles.lh_weights, 1);
    rh_corr = sum(rh_corr(handles.rh_vno) .* handles.rh_weights, 1);
    
    lh_corr(lh_corr == 0) = handles.lh_curv(lh_corr == 0);
    rh_corr(rh_corr == 0) = handles.rh_curv(rh_corr == 0);
    lh_corr(handles.lh_boundary_vec) = 2;
    rh_corr(handles.rh_boundary_vec) = 2;
    lh_corr(1) = 0;
    rh_corr(1) = 0;
else
    lh_corr(lh_corr == 0) = handles.lh_curv(lh_corr == 0);
    rh_corr(rh_corr == 0) = handles.rh_curv(rh_corr == 0);
    lh_corr(1) = 0;
    rh_corr(1) = 0;
end

subplot(2, 2, 1); TrisurfMeshData(handles.lh_avg_mesh, lh_corr); shading interp; colormap(handles.corr_colormap); colorbar
subplot(2, 2, 2); TrisurfMeshData(handles.rh_avg_mesh, rh_corr); shading interp; colormap(handles.corr_colormap); colorbar
subplot(2, 2, 3); TrisurfMeshData(handles.lh_avg_mesh, lh_corr); shading interp; colormap(handles.corr_colormap); colorbar
subplot(2, 2, 4); TrisurfMeshData(handles.rh_avg_mesh, rh_corr); shading interp; colormap(handles.corr_colormap); colorbar

subplot(2, 2, 1); view(-90, 0); colorbar off; colorbar; axis off;
subplot(2, 2, 2); view(90, 0); colorbar off; colorbar; axis off; 
subplot(2, 2, 3); view(90, 0); colorbar off; colorbar; axis off; 
subplot(2, 2, 4); view(-90, 0); colorbar off; colorbar; axis off; 
fprintf(' DONE! \n')


function FS_DistToCortex_Callback(hObject, eventdata, handles)
% hObject    handle to FS_DistToCortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FS_DistToCortex as text
%        str2double(get(hObject,'String')) returns contents of FS_DistToCortex as a double


% --- Executes during object creation, after setting all properties.
function FS_DistToCortex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FS_DistToCortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MNI_DIst2Cortex_Callback(hObject, eventdata, handles)
% hObject    handle to MNI_DIst2Cortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MNI_DIst2Cortex as text
%        str2double(get(hObject,'String')) returns contents of MNI_DIst2Cortex as a double


% --- Executes during object creation, after setting all properties.
function MNI_DIst2Cortex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MNI_DIst2Cortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function vertex = GetVertexFromVoxel(vox, vertex_map)

vertex = vertex_map.vol(vox(2)+1, vox(1)+1, vox(3)+1);

function vertex = GetVertexFromRAS(ras, vertex_map)

vox = CBIG_ConvertRas2Vox([ras(1); ras(2); ras(3)], vertex_map.vox2ras);
vertex = GetVertexFromVoxel(vox, vertex_map);



function rh_vertex_Callback(hObject, eventdata, handles)
% hObject    handle to rh_vertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rh_vertex as text
%        str2double(get(hObject,'String')) returns contents of rh_vertex as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

hemi = 'rh';
vertex = str2double(get(handles.rh_vertex,'String'));

% Update freesurfer vox, ras
set(handles.FS_vox, 'String', [num2str(handles.rh_FS_vox(1, vertex), '%.1f') ', ' num2str(handles.rh_FS_vox(2, vertex), '%.1f') ', ' num2str(handles.rh_FS_vox(3, vertex), '%.1f')]);
set(handles.FS_ras, 'String', [num2str(handles.rh_FS_ras(1, vertex), '%.1f') ', ' num2str(handles.rh_FS_ras(2, vertex), '%.1f') ', ' num2str(handles.rh_FS_ras(3, vertex), '%.1f')]);

% Update FSL_MNI152 vox, ras
set(handles.FSL_MNI152_vox, 'String', [num2str(handles.rh_FSL_MNI152_vox(1, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_vox(2, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_vox(3, vertex), '%.1f')]);
set(handles.FSL_MNI152_ras, 'String', [num2str(handles.rh_FSL_MNI152_ras(1, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_ras(2, vertex), '%.1f') ', ' num2str(handles.rh_FSL_MNI152_ras(3, vertex), '%.1f')]);
vertex = vertex + 10242;

% update figure
handles.vertex = vertex;
guidata(hObject, handles);
UpdateFigures(handles, vertex);

% Update distance.
set(handles.FS_DistToCortex, 'String', 'Dist to cortex : surf mode');
set(handles.FSL_MNI152_DistToCortex, 'String', 'Dist to cortex : surf mode');
set(handles.lh_vertex, 'String', 'lh vertex');


% --- Executes during object creation, after setting all properties.
function rh_vertex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rh_vertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lh_vertex_Callback(hObject, eventdata, handles)
% hObject    handle to lh_vertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lh_vertex as text
%        str2double(get(hObject,'String')) returns contents of lh_vertex as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

hemi = 'lh';
vertex = str2double(get(handles.lh_vertex,'String'));

% Update freesurfer vox, ras
set(handles.FS_vox, 'String', [num2str(handles.lh_FS_vox(1, vertex), '%.1f') ', ' num2str(handles.lh_FS_vox(2, vertex), '%.1f') ', ' num2str(handles.lh_FS_vox(3, vertex), '%.1f')]);
set(handles.FS_ras, 'String', [num2str(handles.lh_FS_ras(1, vertex), '%.1f') ', ' num2str(handles.lh_FS_ras(2, vertex), '%.1f') ', ' num2str(handles.lh_FS_ras(3, vertex), '%.1f')]);

% Update FSL_MNI152 vox, ras
set(handles.FSL_MNI152_vox, 'String', [num2str(handles.lh_FSL_MNI152_vox(1, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_vox(2, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_vox(3, vertex), '%.1f')]);
set(handles.FSL_MNI152_ras, 'String', [num2str(handles.lh_FSL_MNI152_ras(1, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_ras(2, vertex), '%.1f') ', ' num2str(handles.lh_FSL_MNI152_ras(3, vertex), '%.1f')]);

% update figure
handles.vertex = vertex;
guidata(hObject, handles);
UpdateFigures(handles, vertex);

% Update distance.
set(handles.FS_DistToCortex, 'String', 'Dist to cortex : surf mode');
set(handles.FSL_MNI152_DistToCortex, 'String', 'Dist to cortex : surf mode');
set(handles.rh_vertex, 'String', 'lh vertex');

% --- Executes during object creation, after setting all properties.
function lh_vertex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lh_vertex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


