function varargout = Vol2SurfGui(varargin)
% VOL2SURFGUI M-file for Vol2SurfGui.fig
%      VOL2SURFGUI, by itself, creates a new VOL2SURFGUI or raises the existing
%      singleton*.
%
%      H = VOL2SURFGUI returns the handle to a new VOL2SURFGUI or the handle to
%      the existing singleton*.
%
%      VOL2SURFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOL2SURFGUI.M with the given input arguments.
%
%      VOL2SURFGUI('Property','Value',...) creates a new VOL2SURFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Vol2SurfGui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Vol2SurfGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Vol2SurfGui

% Last Modified by GUIDE v2.5 26-Sep-2010 18:04:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Vol2SurfGui_OpeningFcn, ...
                   'gui_OutputFcn',  @Vol2SurfGui_OutputFcn, ...
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


% --- Executes just before Vol2SurfGui is made visible.
function Vol2SurfGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Vol2SurfGui (see VARARGIN)

% Choose default command line output for Vol2SurfGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Vol2SurfGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

% add optional data dir input
[a, hostname] = system('hostname -d');
data_dir = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'vol2surf_gui_data', 'PrecomputedVol2SurfCorrelation', 'LooseCerebellum.dist1.Mask.fwhm6.smooth6.regress1.500sub.1');
p.addParamValue('data_dir', data_dir);

% add option for opengl or zbuffer renderer 
p.addParamValue('renderer', 'zbuffer', @(x)strcmp(x,'opengl') || ...
              strcmp(x,'zbuffer'));

% add optional number of networks. 0 ==> don't load networks
p.addParamValue('num_networks', 0);

% add dataset_num used to estimate the networks.
p.addParamValue('dataset_num', 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provide help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(sum(strcmp(varargin, 'help')) > 0)
   disp(['Vol2SurfGui <data_dir data_dir_string> <renderer renderer_string> <num_networks num_networks> <dataset_num dataset_num>']);
   disp(' ');
    disp('Example Calls: ');
    disp(['    Vol2SurfGui renderer opengl']);
    disp(['    Vol2SurfGui num_networks 7 data_dir ' data_dir]);
    disp(['    Vol2SurfGui renderer zbuffer dataset_num 3 num_networks 7 data_dir ' data_dir]);
    disp(' ');
    disp(['default data_dir = ' data_dir]);
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
data_dir = p.Results.data_dir;
num_networks = p.Results.num_networks;
dataset_num = p.Results.dataset_num;

if(ischar(num_networks))
    num_networks = str2num(num_networks);
end

if(ischar(dataset_num))
    dataset_num = str2num(dataset_num);
end

handles.data_dir = data_dir;
handles.num_networks = num_networks;
handles.dataset_num = dataset_num;
handles.renderer = p.Results.renderer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read mesh
disp('Reading meshes ...');
if(num_networks > 0)
    handles.lh_avg_mesh5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'sphere', 'cortex');
    handles.rh_avg_mesh5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'sphere', 'cortex'); 
    handles.lh_avg_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
    handles.rh_avg_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex'); 
    handles.lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'inflated', 'cortex');
    handles.rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'inflated', 'cortex'); 
%    handles.lh_avg_mesh_inflated5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
%    handles.rh_avg_mesh_inflated5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex'); 
else
    handles.lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
    handles.rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex'); 
    handles.lh_avg_mesh5 = handles.lh_avg_mesh;
    handles.rh_avg_mesh5 = handles.rh_avg_mesh;
%    handles.lh_avg_mesh_inflated5 = handles.lh_avg_mesh;
%    handles.rh_avg_mesh_inflated5 = handles.rh_avg_mesh;
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


basedir = basename(data_dir);
% read template
template_path = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', 'FS_nonlinear_volumetric_space_4.5');
MNI_path = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'volume', 'FSL_MNI152_FS4.5.0');

template_file = fullfile(template_path, 'Brain.1.100.1.download_subjects_sorted1002.txt.2mm.nii.gz');
disp(['Reading template: ' template_file]);
handles.template = MRIread(template_file);
handles.template.vol = handles.template.vol/max(handles.template.vol(:));

MNI_file = fullfile(MNI_path, 'mri', 'norm.mgz');
handles.FSL_MNI152_template = MRIread(MNI_file);

% read mask
mask_str_index = strfind(basedir, 'Mask.');
mask_file = fullfile(template_path, [basedir(1:mask_str_index+4) 'GCA.t0.5.nii.gz']);
disp(['Reading mask: ' mask_file]);
handles.mask = MRIread(mask_file);

% read parcellations
if(num_networks > 0)
    
    % grab key words from
    keys{1} = ['.' strtok(basedir, '.') '.'];
    remainder = basedir;
    count = 1;
    while(~isempty(remainder))
        [cand_keys{count}, remainder] = strtok(remainder, '.');
        count = count + 1;
    end

    count = 1;
    possible_keys = {'dist', 'fwhm', 'smooth', 'regress'};
    for i = 1:length(possible_keys)
        for j = 1:length(cand_keys)
            key_index = strfind(cand_keys{j}, possible_keys{i});
            if(~isempty(key_index))
                count = count + 1;
                keys{count} = ['.' cand_keys{j} '.'];
            end
        end
    end

    % Reading volume segmentation
    num_networks_str = num2str(num_networks, '%03d');
    subcort_path = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'vol2surf_gui_data', 'subcortical_cluster_results', ['cluster' num_networks_str]);

    curr_pwd = pwd;
    cd(subcort_path);
    if(dataset_num == 2)
        cand = dir('*.500sub.2.*.parcel.nii.gz');
    elseif(dataset_num == 1)
        cand = dir('*.500sub.1.*.parcel.nii.gz');
    elseif(dataset_num == 3)
        cand = dir('*.1000sub.*.parcel.nii.gz');
    else
        error(['Does not handle data set num = ' num2str(dataset_num)]);
    end
    cd(curr_pwd);

    index = [];
    for i = 1:length(cand)
        cand_bool = 1;
        for j = 1:length(keys)
            cand_bool = cand_bool & ~isempty(strfind(cand(i).name, keys{j}));
        end
        if(cand_bool)
            index = [index; i];
        end
    end

    if(length(index) > 1)
        warning('There are more than 1 subcort segmentation file satisfying criteria. Picking first one');
    elseif(isempty(index))
        error('There are no segmentation file satisfying criteria.');
    end
    subcort_network_file = fullfile(subcort_path, cand(index(1)).name);
    disp(['Reading volume segmentation: ' subcort_network_file]);
    handles.subcort_network = MRIread(subcort_network_file);
    handles.subcort_network.vol = handles.subcort_network.vol + 1; % so network ranges from 1 to network+1 
    % background replaced with anatomy, also networks should theoretically range from 2 to network+1
    handles.subcort_network.vol(handles.mask.vol == 0) = handles.template.vol(handles.mask.vol == 0);
    
    % Reading volume silhouette
    curr_pwd = pwd;
    cd(subcort_path);
    if(dataset_num == 2)
        cand = dir('*.500sub.2.*.cort_silhouette.nii.gz');
    elseif(dataset_num == 1)
        cand = dir('*.500sub.1.*.cort_silhouette.nii.gz');
    elseif(dataset_num == 3)
        cand = dir('*.1000sub.*.cort_silhouette.nii.gz');
    else
        error(['Does not handle data set num = ' num2str(dataset_num)]);
    end
    cd(curr_pwd);

    index = [];
    for i = 1:length(cand)
        cand_bool = 1;
        for j = 1:length(keys)
            cand_bool = cand_bool & ~isempty(strfind(cand(i).name, keys{j}));
        end
        if(cand_bool)
            index = [index; i];
        end
    end

    if(length(index) > 1)
        warning('There are more than 1 subcort silhouette file satisfying criteria. Picking first one');
    elseif(isempty(index))
        error('There are no silhouette file satisfying criteria');
    end
    subcort_sil_file = fullfile(subcort_path, cand(index(1)).name);
    disp(['Reading subcortical silhouette: ' subcort_sil_file]);
    handles.subcort_sil = MRIread(subcort_sil_file);
    handles.subcort_sil.vol = handles.subcort_sil.vol + 2; % so sil ranges from 1 to 3
    handles.subcort_sil.vol(handles.mask.vol == 0) = handles.template.vol(handles.mask.vol == 0); %background replaced with anatomy
    
    % Reading surface parcellations
    num_networks_str_prev = num2str(num_networks-1, '%03d');
    data_path = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'extras', 'surf2surf_gui_data', 'cluster_results', ['cluster' num_networks_str]);
    colortable_path = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'utilities');

    if(dataset_num == 2)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.500sub.2.ref500sub.1.mat']);
    elseif(dataset_num == 1)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '_ref' num_networks_str_prev '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.500sub.1.mat']);
    elseif(dataset_num == 3)
        network_file = fullfile(data_path, ['Cluster' num_networks_str '.s00.tries1000.rand000.znorm0.profile.fsaverage3.smooth6.s5.t0.1.1000sub.ref500sub.1.mat']);
    else
        error(['Does not handle data set num = ' num2str(dataset_num)]);
    end
    disp(['Reading surface parcellations: ' network_file]);
    load(network_file);
    handles.lh_labels = lh_labels';
    handles.rh_labels = rh_labels';
    
    % upsample labels to high resolution
    handles.lh_labels7 = MARS_NNInterpolate(handles.lh_avg_mesh7.vertices, handles.lh_avg_mesh5, handles.lh_labels);
    handles.rh_labels7 = MARS_NNInterpolate(handles.rh_avg_mesh7.vertices, handles.rh_avg_mesh5, handles.rh_labels);

    load(fullfile(colortable_path, 'freesurfer_color.mat'));
    handles.freesurfer_color = freesurfer_color(1:num_networks+1, :);

    % precompute interpolation weights
    [tmp, tmp, NF] = MARS_linearInterpolate(handles.lh_avg_mesh7.vertices, handles.lh_avg_mesh5, handles.lh_avg_mesh5.data(1, :));
    [handles.lh_weights, handles.lh_vno] = MARS_computeWeightsForInterpolation(handles.lh_avg_mesh5, handles.lh_avg_mesh7.vertices, NF);
    [tmp, tmp, NF] = MARS_linearInterpolate(handles.rh_avg_mesh7.vertices, handles.rh_avg_mesh5, handles.rh_avg_mesh5.data(1, :));
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
    
    % create colormap for parcellations
    network_colormap = gray(128);
    num_colors = size(network_colormap, 1);
    network_colormap = [network_colormap; repmat(network_colormap(end, :), num_colors/2, 1)];
    for i = 2:size(handles.freesurfer_color, 1)-1
        network_colormap = [network_colormap; repmat(handles.freesurfer_color(i, :), num_colors, 1)];
    end
    handles.network_colormap = [network_colormap; repmat(handles.freesurfer_color(end, :), num_colors/2, 1)];

    % create colormap for silhouette
    max_sil = max(handles.subcort_sil.vol(handles.mask.vol == 1));
    min_sil = min(handles.subcort_sil.vol(handles.mask.vol == 1));
    jet_colormap = jet(128);
    num_high_colors = round(size(jet_colormap, 1) / (max_sil - min_sil) * (3 - max_sil));
    num_low_colors  = round(size(jet_colormap, 1) / (max_sil - min_sil) * (min_sil - 1)/2); %only half
    padded_jet_colormap = [repmat(jet_colormap(1, :), num_low_colors, 1) ; jet_colormap; repmat(jet_colormap(end, :), num_high_colors, 1)];
    
    total_gray_colors = round(size(jet_colormap, 1) / (max_sil - min_sil));
    num_remain_colors = total_gray_colors * 2 - size(padded_jet_colormap, 1);
    tmp = gray(total_gray_colors);
    handles.sil_colormap = [tmp; repmat(tmp(end, :), num_remain_colors, 1); padded_jet_colormap];
end
% create new colormap;
handles.corr_colormap = CreateCorrColormap(handles.curv_gray_scale, low_threshold, high_threshold);

% dim the template that is not in the mask.
handles.template.vol(handles.mask.vol == 0) = handles.template.vol(handles.mask.vol == 0)*0.25;

% reading coordinate systems
disp('Loading Coordinate Transforms');
coord_path = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'transforms', 'CorrespondenceFreeSurferVolSurfSpace_Buckner2011');
load(fullfile(coord_path, 'coord_MNI152_FSvol', 'FS2mm2MNI1mm.vox.mat')); handles.FS2FSL_MNI152_vox = vox;
load(fullfile(coord_path, 'coord_MNI152_FSvol', 'FS2mm2MNI1mm.ras.mat')); handles.FS2FSL_MNI152_ras = ras;
load(fullfile(coord_path, 'coord_MNI152_FSvol', 'MNI1mm2FS2mm.vox.mat')); handles.FSL_MNI1522FS_vox = vox;
load(fullfile(coord_path, 'coord_MNI152_FSvol', 'MNI1mm2FS2mm.ras.mat')); handles.FSL_MNI1522FS_ras = ras;


% create mapping from original to new volume
im = handles.template.vol(:, :, 1);
handles.orig_coronal_im = reshape(1:numel(im), size(im, 1), size(im, 2));
handles.new_coronal_im = handles.orig_coronal_im(:, end:-1:1);

im = squeeze(handles.template.vol(1, :, :));
handles.orig_axial_im = reshape(1:numel(im), size(im, 1), size(im, 2));
handles.new_axial_im = handles.orig_axial_im';
handles.new_axial_im = handles.new_axial_im(end:-1:1, end:-1:1);

im = squeeze(handles.template.vol(:, 1, :));
handles.orig_sagittal_im = reshape(1:numel(im), size(im, 1), size(im, 2));
handles.new_sagittal_im = handles.orig_axial_im;

% Compute best slice to display.
[Y, coronal_slice] = max(sum(squeeze(sum(handles.mask.vol, 1)), 1));
[Y, axial_slice] = max(sum(squeeze(sum(handles.mask.vol, 2)), 2));
[Y, sagittal_slice] = max(sum(squeeze(sum(handles.mask.vol, 1)), 2));
FS_vox = [sagittal_slice-1; axial_slice-1 ; coronal_slice-1];
DrawSlices(handles, FS_vox);

% set coordinates strings
UpdateCoordinatesFromFSvox(FS_vox, hObject, handles)

% Set up data cursor structure
dcm_obj = datacursormode(handles.figure1);
handles.dcm_obj = dcm_obj;
set(dcm_obj,'UpdateFcn',@myupdatefcn)

% initialize voxels
handles.voxel = [];
handles.curr_voxel = FS_vox;



global super_handles; super_handles = handles;
guidata(hObject, handles);





function txt = myupdatefcn(empt, event_obj)
% Customizes text of data tips

global super_handles;
handles = guidata(super_handles.figure1);

if(~isempty(find(findall(super_handles.coronal) == get(event_obj, 'Target'), 1)))
    slice_type = 'coronal';
    slice = handles.curr_voxel(3) + 1;
    orig_im = handles.orig_coronal_im;
    new_im = handles.new_coronal_im;
elseif(~isempty(find(findall(super_handles.axial) == get(event_obj, 'Target'), 1)))
    slice_type = 'axial';
    slice = handles.curr_voxel(2) + 1;
    orig_im = handles.orig_axial_im;
    new_im = handles.new_axial_im;
elseif(~isempty(find(findall(super_handles.sagittal) == get(event_obj, 'Target'), 1)))    
    slice_type = 'sagittal';
    slice = handles.curr_voxel(1) + 1;
    orig_im = handles.orig_sagittal_im;
    new_im = handles.new_sagittal_im;
else
    warning('Please select point from coronal, axial or sagittal slices');
end

[x, y, z] = FindVoxelCoordinates(slice_type, slice, round(get(event_obj,'Position')'), orig_im, new_im);
FS_ras = CBIG_ConvertVox2Ras([x;y;z], handles.template.vox2ras);

if(handles.num_networks > 0)
    network = handles.subcort_network.vol(y+1,x+1,z+1);
    if(network > 1)
        txt = {['FS_X: ',num2str(FS_ras(1))],...
            ['FS_Y: ', num2str(FS_ras(2))], ...
            ['FS_Z: ', num2str(FS_ras(3))], ...
            ['network: ', num2str(network - 1)], ...
            ['sil: ', num2str(handles.subcort_sil.vol(y+1,x+1,z+1) - 2)]};
    else
        txt = {['FS_X: ',num2str(FS_ras(1))],...
        ['FS_Y: ', num2str(FS_ras(2))], ...
        ['FS_Z: ', num2str(FS_ras(3))]};
    end
else
    txt = {['FS_X: ',num2str(FS_ras(1))],...
        ['FS_Y: ', num2str(FS_ras(2))], ...
        ['FS_Z: ', num2str(FS_ras(3))]};
end

function [x, y, z] = FindVoxelCoordinates(slice_type, slice, Position, orig_im, new_im)

% note: this is NOT ras
[orig_row, orig_col] = find(orig_im == new_im(Position(2), Position(1)));

if(strcmp(slice_type, 'coronal'))
    x = orig_col-1;
    y = orig_row-1;
    z =  slice-1;
elseif(strcmp(slice_type, 'axial'))
    x = orig_row-1;
    y = slice-1;
    z = orig_col-1 ;
elseif(strcmp(slice_type, 'sagittal'))
    x = slice-1;
    y = orig_row-1;
    z = orig_col-1;
end


function DrawSlices(handles, FS_vox)

%handles.curr_voxel
if(get(handles.Anatomy_mode, 'Value'))
    image_type = 'template';
elseif(get(handles.Network_mode, 'Value'))
    image_type = 'network';
elseif(get(handles.Confidence_mode, 'Value'))
    image_type = 'sil';
end

DrawCoronalSlice(handles, image_type, FS_vox(3)+1);
DrawAxialSlice(handles, image_type, FS_vox(2)+1);
DrawSagittalSlice(handles, image_type, FS_vox(1)+1);




function DrawCoronalSlice(handles, image_type, slice)

XLim = [25.5 103.5];
YLim = [20 110];

axes(handles.coronal);
if(strcmpi(image_type, 'template'))
    im = squeeze(handles.template.vol(:, :, slice));
    im = im(:, end:-1:1); % transform to neurological LR
    
    im(1, 1, 1) = 0; % ensure we hit the entire range
    im(end, end, end) = 1; %ensure we hit entire range
    imagesc(im, [0 1]);
    colormap(gray(64))
    
elseif(strcmpi(image_type, 'network'))
    im = squeeze(handles.subcort_network.vol(:, :, slice));
    im = im(:, end:-1:1); % transform to neurological LR
    
    im(1, 1, 1) = handles.num_networks+1; % ensure we hit the entire range
    im(end, end, end) = 0; %ensure we hit entire range
    imagesc(im, [0 handles.num_networks+1]);
    colormap(handles.network_colormap);
    
elseif(strcmpi(image_type, 'sil'))
    im = squeeze(handles.subcort_sil.vol(:, :, slice));
    im = im(:, end:-1:1); % transform to neurological LR
    
    im(1, 1, 1) = 3; % ensure we hit the entire range
    im(end, end, end) = 0; %ensure we hit entire range
    imagesc(im, [0 3]);
    colormap(handles.sil_colormap);
else
    error(['Does not handle image type: ' image_type]);
end    

set(handles.coronal, 'XLim', XLim);
set(handles.coronal, 'YLim', YLim);
axis off;



function DrawAxialSlice(handles, image_type, FS_slice)

XLim = [25.5 103.5];
YLim = [25 121];

axes(handles.axial);
if(strcmpi(image_type, 'template'))
    im = squeeze(handles.template.vol(FS_slice, :, :));
    im = im';
    im = im(end:-1:1, end:-1:1);
    
    im(1, 1, 1) = 0; % ensure we hit the entire range
    im(end, end, end) = 1; %ensure we hit entire range
    imagesc(im, [0 1]);
    colormap(gray(64))
    
elseif(strcmpi(image_type, 'network'))
    im = squeeze(handles.subcort_network.vol(FS_slice, :, :));
    im = im';
    im = im(end:-1:1, end:-1:1);
    
    im(1, 1, 1) = handles.num_networks+1; % ensure we hit the entire range
    im(end, end, end) = 0; %ensure we hit entire range
    imagesc(im, [0 handles.num_networks+1]);
    colormap(handles.network_colormap);
    
elseif(strcmpi(image_type, 'sil'))
    im = squeeze(handles.subcort_sil.vol(FS_slice, :, :));
    im = im';
    im = im(end:-1:1, end:-1:1);
    
    im(1, 1, 1) = 3; % ensure we hit the entire range
    im(end, end, end) = 0; %ensure we hit entire range
    imagesc(im, [0 3]);
    colormap(handles.sil_colormap);
else
    error(['Does not handle image type: ' image_type]);
end    

set(handles.axial, 'XLim', XLim);
set(handles.axial, 'YLim', YLim);
axis off;

function DrawSagittalSlice(handles, image_type, FS_slice)

XLim = [8 104];
YLim = [20 110];

axes(handles.sagittal);
if(strcmpi(image_type, 'template'))
    im = squeeze(handles.template.vol(:, FS_slice, :));
    
    im(1, 1, 1) = 0; % ensure we hit the entire range
    im(end, end, end) = 1; %ensure we hit entire range
    imagesc(im, [0 1]);
    colormap(gray(64))
    
elseif(strcmpi(image_type, 'network'))
    im = squeeze(handles.subcort_network.vol(:, FS_slice, :));
    
    im(1, 1, 1) = handles.num_networks+1; % ensure we hit the entire range
    im(end, end, end) = 0; %ensure we hit entire range
    imagesc(im, [0 handles.num_networks+1]);
    colormap(handles.network_colormap);
    
elseif(strcmpi(image_type, 'sil'))
    im = squeeze(handles.subcort_sil.vol(:, FS_slice, :));
    
    im(1, 1, 1) = 3; % ensure we hit the entire range
    im(end, end, end) = 0; %ensure we hit entire range
    imagesc(im, [0 3]);
    colormap(handles.sil_colormap);
else
    error(['Does not handle image type: ' image_type]);
end    

set(handles.sagittal, 'XLim', XLim);
set(handles.sagittal, 'YLim', YLim);
axis off;



% --- Outputs from this function are returned to the command line.
function varargout = Vol2SurfGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

% Get default command line output from handles structure
varargout{1} = handles.output;

if(handles.help)
    delete(handles.figure1);
end

function FS_vox_Callback(hObject, eventdata, handles)
% hObject    handle to FS_vox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

% Hints: get(hObject,'String') returns contents of FS_vox as text
%        str2double(get(hObject,'String')) returns contents of FS_vox as a double

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

vox = round(str2num(get(hObject, 'String')));
handles.curr_voxel = [vox(1); vox(2); vox(3)];

UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);

DrawSlices(handles, handles.curr_voxel);

if(get(handles.SelectVoxel, 'Value'))
    success = UpdateFigures(handles, vox);

    if(success)
        handles = guidata(gcbo);
        handles.voxel = handles.curr_voxel;
        guidata(hObject, handles);
        
        sil = handles.subcort_sil.vol(handles.curr_voxel(2)+1,handles.curr_voxel(1)+1,handles.curr_voxel(3)+1) - 2
    end
end

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
% handles    structure with handles and user data (see guidata)

% Hints: get(hObject,'String') returns contents of FS_ras as text
%        str2double(get(hObject,'String')) returns contents of FS_ras as a double
delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

% this is necessary because we are in 2mm space, so we want to round off to
% the nearest even number
ras = round(str2num(get(hObject, 'String'))/2)*2;
ras = [ras(1); ras(2); ras(3)];
handles.curr_voxel = CBIG_ConvertRas2Vox(ras, handles.template.vox2ras);

UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);

DrawSlices(handles, handles.curr_voxel);

if(get(handles.SelectVoxel, 'Value'))
    success = UpdateFigures(handles, handles.curr_voxel);

    if(success)
        handles = guidata(gcbo);
        handles.voxel = handles.curr_voxel;
        guidata(hObject, handles);
        
        sil = handles.subcort_sil.vol(handles.curr_voxel(2)+1,handles.curr_voxel(1)+1,handles.curr_voxel(3)+1) - 2
    end
end

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
% handles    structure with handles and user data (see guidata)

% Hints: get(hObject,'String') returns contents of FSL_MNI152_vox as text
%        str2double(get(hObject,'String')) returns contents of FSL_MNI152_vox as a double
delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

vox = round(str2num(get(hObject, 'String')));
handles.curr_voxel = GetFSVoxFromMNIVox([vox(1); vox(2); vox(3)], handles.FSL_MNI1522FS_vox);

UpdateCoordinatesFromFSL_MNI152_vox([vox(1); vox(2); vox(3)], hObject, handles);
handles = guidata(gcbo);
guidata(hObject, handles);

DrawSlices(handles, handles.curr_voxel);

if(get(handles.SelectVoxel, 'Value'))
    success = UpdateFigures(handles, handles.curr_voxel);

    if(success)
        handles = guidata(gcbo);
        handles.voxel = handles.curr_voxel;
        guidata(hObject, handles);
        
        sil = handles.subcort_sil.vol(handles.curr_voxel(2)+1,handles.curr_voxel(1)+1,handles.curr_voxel(3)+1) - 2
    end
end

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
% handles    structure with handles and user data (see guidata)

% Hints: get(hObject,'String') returns contents of FSL_MNI152_ras as text
%        str2double(get(hObject,'String')) returns contents of FSL_MNI152_ras as a double

ras = round(str2num(get(hObject, 'String')));
vox = CBIG_ConvertRas2Vox([ras(1); ras(2); ras(3)], handles.FSL_MNI152_template.vox2ras);
handles.curr_voxel = GetFSVoxFromMNIVox(vox, handles.FSL_MNI1522FS_vox);

UpdateCoordinatesFromFSL_MNI152_vox(vox, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);

DrawSlices(handles, handles.curr_voxel);

if(get(handles.SelectVoxel, 'Value'))
    success = UpdateFigures(handles, handles.curr_voxel);

    if(success)
        handles = guidata(gcbo);
        handles.voxel = handles.curr_voxel;
        guidata(hObject, handles);
        
        sil = handles.subcort_sil.vol(handles.curr_voxel(2)+1,handles.curr_voxel(1)+1,handles.curr_voxel(3)+1) - 2
    end
end

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
% handles    structure with handles and user data (see guidata)

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
if(~isempty(handles.voxel))
    UpdateFigures(handles, handles.voxel);
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
% handles    structure with handles and user data (see guidata)

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
if(~isempty(handles.voxel))
    UpdateFigures(handles, handles.voxel);
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


% --- Executes on button press in Posterior.
function Posterior_Callback(hObject, eventdata, handles)
% hObject    handle to Posterior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
set(handles.Navigate, 'Value', 1);
set(handles.SelectVoxel, 'Value', 0);


handles.curr_voxel(3) = max(handles.curr_voxel(3) - 1, 0);
%handles.curr_voxel
DrawSlices(handles, handles.curr_voxel);

UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);


% --- Executes on button press in Anterior.
function Anterior_Callback(hObject, eventdata, handles)
% hObject    handle to Anterior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
set(handles.Navigate, 'Value', 1);
set(handles.SelectVoxel, 'Value', 0);


handles.curr_voxel(3) = min(handles.curr_voxel(3) + 1, 127);
%handles.curr_voxel
DrawSlices(handles, handles.curr_voxel);
    
UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);


% --- Executes on button press in Superior.
function Superior_Callback(hObject, eventdata, handles)
% hObject    handle to Superior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
set(handles.Navigate, 'Value', 1);
set(handles.SelectVoxel, 'Value', 0);


handles.curr_voxel(2) = max(handles.curr_voxel(2) - 1, 0);
%handles.curr_voxel
DrawSlices(handles, handles.curr_voxel);
    
UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);

% --- Executes on button press in Inferior.
function Inferior_Callback(hObject, eventdata, handles)
% hObject    handle to Inferior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
set(handles.Navigate, 'Value', 1);
set(handles.SelectVoxel, 'Value', 0);


handles.curr_voxel(2) = min(handles.curr_voxel(2) + 1, 127);
%handles.curr_voxel
DrawSlices(handles, handles.curr_voxel);
    
UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);

% --- Executes on button press in Left.
function Left_Callback(hObject, eventdata, handles)
% hObject    handle to Left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
set(handles.Navigate, 'Value', 1);
set(handles.SelectVoxel, 'Value', 0);


handles.curr_voxel(1) = min(handles.curr_voxel(1) + 1, 127);
%handles.curr_voxel
DrawSlices(handles, handles.curr_voxel);
    
UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles)
handles = guidata(gcbo);
guidata(hObject, handles);

% --- Executes on button press in Right.
function Right_Callback(hObject, eventdata, handles)
% hObject    handle to Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));
set(handles.Navigate, 'Value', 1);
set(handles.SelectVoxel, 'Value', 0);


handles.curr_voxel(1) = max(handles.curr_voxel(1) - 1, 0);
%handles.curr_voxel
DrawSlices(handles, handles.curr_voxel);
    
UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles);
handles = guidata(gcbo);
guidata(hObject, handles);






function MNI_vox = GetMNIVoxFromFSVox(vox, FS2FSL_MNI152_vox)

index = sub2ind([128 128 128], vox(2)+1, vox(1)+1, vox(3) + 1);
MNI_vox = FS2FSL_MNI152_vox(:, index);

function FS_vox = GetFSVoxFromMNIVox(vox, FSL_MNI1522FS_vox)

index = sub2ind([256 256 256], vox(2)+1, vox(1)+1, vox(3) + 1);
FS_vox = FSL_MNI1522FS_vox(:, index);

function UpdateCoordinatesFromFSvox(FS_vox, hObject, handles)

FS_ras = CBIG_ConvertVox2Ras(FS_vox, handles.template.vox2ras);
FSL_MNI152_vox = GetMNIVoxFromFSVox(FS_vox, handles.FS2FSL_MNI152_vox);
FSL_MNI152_ras = CBIG_ConvertVox2Ras(FSL_MNI152_vox, handles.FSL_MNI152_template.vox2ras);

set(handles.FS_vox, 'String', [num2str(FS_vox(1)) ', ' num2str(FS_vox(2)) ', ' num2str(FS_vox(3))]);
set(handles.FS_ras, 'String', [num2str(FS_ras(1)) ', ' num2str(FS_ras(2)) ', ' num2str(FS_ras(3))]);
set(handles.FSL_MNI152_vox, 'String', [num2str(FSL_MNI152_vox(1)) ', ' num2str(FSL_MNI152_vox(2)) ', ' num2str(FSL_MNI152_vox(3))]);
set(handles.FSL_MNI152_ras, 'String', [num2str(FSL_MNI152_ras(1)) ', ' num2str(FSL_MNI152_ras(2)) ', ' num2str(FSL_MNI152_ras(3))]);

guidata(hObject, handles);

function UpdateCoordinatesFromFSL_MNI152_vox(FSL_MNI152_vox, hObject, handles)

FSL_MNI152_ras = CBIG_ConvertVox2Ras(FSL_MNI152_vox, handles.FSL_MNI152_template.vox2ras);
FS_vox = GetFSVoxFromMNIVox(FSL_MNI152_vox, handles.FSL_MNI1522FS_vox);
FS_ras = CBIG_ConvertVox2Ras(FS_vox, handles.template.vox2ras);

set(handles.FS_vox, 'String', [num2str(FS_vox(1)) ', ' num2str(FS_vox(2)) ', ' num2str(FS_vox(3))]);
set(handles.FS_ras, 'String', [num2str(FS_ras(1)) ', ' num2str(FS_ras(2)) ', ' num2str(FS_ras(3))]);
set(handles.FSL_MNI152_vox, 'String', [num2str(FSL_MNI152_vox(1)) ', ' num2str(FSL_MNI152_vox(2)) ', ' num2str(FSL_MNI152_vox(3))]);
set(handles.FSL_MNI152_ras, 'String', [num2str(FSL_MNI152_ras(1)) ', ' num2str(FSL_MNI152_ras(2)) ', ' num2str(FSL_MNI152_ras(3))]);

guidata(hObject, handles);

% --- Executes on button press in Navigate.
function Navigate_Callback(hObject, eventdata, handles)
% hObject    handle to Navigate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

% Hint: get(hObject,'Value') returns toggle state of Navigate
set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','off');

delete(findall(handles.figure1,'Type','hggroup','HandleVisibility','off'));

% --- Executes on button press in SelectVoxel.
function SelectVoxel_Callback(hObject, eventdata, handles)
% hObject    handle to SelectVoxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see guidata)

% Hint: get(hObject,'Value') returns toggle state of SelectVoxel

while(get(hObject, 'Value'))
  set(handles.dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on');
  
  disp('Select voxel and press enter');  
  pause                            % Wait while the user does this.
  cursor_info = getCursorInfo(handles.dcm_obj);
  if(~isempty(cursor_info))
      handles = guidata(gcbo); % this is necessary because stuff like the colormap might have changed between clicking vertices

      if(~isempty(find(findall(handles.coronal) == cursor_info.Target, 1)))
          slice_type = 'coronal';
          slice = handles.curr_voxel(3) + 1;
          orig_im = handles.orig_coronal_im;
          new_im = handles.new_coronal_im;
      elseif(~isempty(find(findall(handles.axial) == cursor_info.Target, 1)))
          slice_type = 'axial';
          slice = handles.curr_voxel(2) + 1;
          orig_im = handles.orig_axial_im;
          new_im = handles.new_axial_im;
      elseif(~isempty(find(findall(handles.sagittal) == cursor_info.Target, 1)))
          slice_type = 'sagittal';
          slice = handles.curr_voxel(1) + 1;
          orig_im = handles.orig_sagittal_im;
          new_im = handles.new_sagittal_im;
      else
          warning('Please select point from coronal, axial or sagittal slices');
      end

      [x, y, z] = FindVoxelCoordinates(slice_type, slice, round(cursor_info.Position'), orig_im, new_im);
      FS_vox = [x;y;z];
      
      success = UpdateFigures(handles, FS_vox);
      
      if(success)
          handles = guidata(gcbo);
          handles.voxel = FS_vox;
          handles.curr_voxel = FS_vox;
          DrawSlices(handles, handles.curr_voxel);
          guidata(hObject, handles);
          
          UpdateCoordinatesFromFSvox(handles.curr_voxel, hObject, handles);
      end
        
  end
end

function success = UpdateFigures(handles, voxel)

h = figure(2); set(gcf, 'renderer', handles.renderer);
pos = get(h, 'Position');
pos(3) = 1000;
pos(4) = 700;
set(h, 'Position', pos);

lh_corr_file = fullfile(handles.data_dir, ['lh.' num2str(voxel(1)) '_' num2str(voxel(2)) '_' num2str(voxel(3)) '.corr.nii.gz']);
disp(lh_corr_file);
if(exist(lh_corr_file, 'file'))
    lh_corr = MRIread(lh_corr_file);
    lh_corr = lh_corr.vol(:);
else
    disp(['Correlation file not found: ' lh_corr_file]);
    success = 0;
    return
end

rh_corr_file = fullfile(handles.data_dir, ['rh.' num2str(voxel(1)) '_' num2str(voxel(2)) '_' num2str(voxel(3)) '.corr.nii.gz']);
disp(rh_corr_file);
if(exist(rh_corr_file, 'file'))
    rh_corr = MRIread(rh_corr_file);
    rh_corr = rh_corr.vol(:);
else
    disp(['Correlation file not found: ' rh_corr_file]);
    success = 0;
    return
end

low_threshold = str2double(get(handles.low_threshold,'String'));
high_threshold = str2double(get(handles.high_threshold,'String'));

lh_corr(handles.lh_avg_mesh5.MARS_label == 1) = 0; lh_corr(lh_corr < low_threshold) = 0; lh_corr(lh_corr > high_threshold) = high_threshold;
rh_corr(handles.rh_avg_mesh5.MARS_label == 1) = 0; rh_corr(rh_corr < low_threshold) = 0; rh_corr(rh_corr > high_threshold) = high_threshold;

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
fprintf(' DONE! \n');
success = 1;

% --- Executes on button press in Anatomy_mode.
function Anatomy_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Anatomy_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Anatomy_mode
if(handles.num_networks > 0)
    DrawSlices(handles, handles.curr_voxel);
end
    
% --- Executes on button press in Network_mode.
function Network_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Network_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Network_mode
if(handles.num_networks > 0)
    DrawSlices(handles, handles.curr_voxel);
else
    set(handles.Network_mode, 'Value', 0);
    set(handles.Anatomy_mode, 'Value', 1);
end


% --- Executes on button press in Confidence_mode.
function Confidence_mode_Callback(hObject, eventdata, handles)
% hObject    handle to Confidence_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Confidence_mode
if(handles.num_networks > 0)
    DrawSlices(handles, handles.curr_voxel);
else
    set(handles.Confidence_mode, 'Value', 0);
    set(handles.Anatomy_mode, 'Value', 1);
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

