function [lh_data_resample, rh_data_resample] = CBIG_resample_fslr(lh_data, rh_data, orig_mesh, resample_mesh, type_of_data, folder_to_write, registration_version)

%[lh_data_resample, rh_data_resample] = CBIG_resample_fslr(lh_data, rh_data, orig_mesh, resample_mesh, type_of_data, folder_to_write, registration_version)
%
% This function resamples label/metric data in
% fs_LR_32k to fs_LR_164k or fs_LR_164k to fs_LR_32k. 
% Input:
%      -lh_data, rh_data:
%       a N x 1 vector, where N is the number of vertices. ?h_data can
%       be metric data (float) or label data (integer). 
%
%      -orig_mesh:
%       'fs_LR_32k' or 'fs_LR_164k'. Mesh name of lh_data, 
%       rh_data.
%
%      -resample_mesh:
%       'fs_LR_32k' or 'fs_LR_164k'. Mesh name of lh_data_resample, 
%       rh_data_resample. The target mesh of the resampling.
%
%      -type_of_data:
%       'metric' or 'label'. If ?h_data is a metric data (float), then
%       type_of_data should be set to 'metric'. IF ?h_data is a label 
%       data (integer), then type_of_data should be set to 'label'.
%
%      -folder_to_write:
%       output path. e.g. '/data/Resample_fsLR'.
%       Note: folder_to_write must be a non-existent directory or folder, 
%             as this folder will be removed when projection is done.
%
%      -registration_version:
%       '20170508' or '20160827'. Human Connectome Project (HCP) group  
%       provides two atlas-to-atlas registration versions. If user
%       doesn't pass in registration_version, default version is 
%       '20170508'.
%
% Output:
%      -lh_data_resample, rh_data_resample:
%       output data in <resample_mesh>.
%
% Example:
% folder_to_write='/data/users/rkong/storage/ruby/data/HCP_relevant/Mapping_FS2fsLR'
% [lh_label_fsLR_164k,rh_label_fsLR_164k]=CBIG_resample_fslr(lh_data_32k,rh_data_32k,'fs_LR_32k','fs_LR_164k','label',folder_to_write,'20170508');

% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(lh_data,2) ~= 1
    error('Input argument ''lh_data'' should be a column vector');
end
if size(rh_data,2) ~= 1
    error('Input argument ''rh_data'' should be a column vector');
end

if(nargin<5)% if you dont set your own write folder
    error('Not enough inputs')
end
if(nargin<6)% if you dont set the atlas-to-atlas registration version
    registration_version = '20170508';
end
if(exist(folder_to_write,'dir') ~= 0) % if the output folder exists
    error('The output folder already exists, please set a non-existent output folder!');
else
    mkdir(folder_to_write);
end


%% Save input ?h_data as a cifti file
% set gifti file extension
if(strcmp(type_of_data, 'label'))
    fsLR_extension = 'label.gii'; % if type_of_data is label, extension of gifti file is .label.gii
elseif(strcmp(type_of_data, 'metric'))% assuming float data
    fsLR_extension = 'func.gii'; % if type_of_data is metric, extension of gifti file is .func.gii
else
    error('unknown type of data')
end

% get cifti func file structure
lh_fsLR_orig = gifti(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', orig_mesh, 'example_func', 'Parcels_L.func.gii')); 
rh_fsLR_orig = gifti(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', orig_mesh, 'example_func', 'Parcels_L.func.gii')); 

lh_fsLR_orig.cdata = lh_data;
rh_fsLR_orig.cdata = rh_data;

% save as a gifti file
save(lh_fsLR_orig, fullfile(folder_to_write, [type_of_data, '_L.', orig_mesh, '.', fsLR_extension]), 'Base64Binary');
save(rh_fsLR_orig, fullfile(folder_to_write, [type_of_data, '_R.', orig_mesh, '.', fsLR_extension]), 'Base64Binary');

%% Perform resampling using wb_command with CBIG_resample_fslr_data.sh

system([fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'fslr_matlab', 'CBIG_resample_fslr_data.sh'), ' ', orig_mesh, ' ', resample_mesh, ' ', folder_to_write, ' ', type_of_data, ' ', registration_version]);

%% Output fs_LR_32k/fs_LR_164k data
lh_data_resample_gifti = gifti(fullfile(folder_to_write, [type_of_data, '_L.', resample_mesh, '.', fsLR_extension]));
rh_data_resample_gifti = gifti(fullfile(folder_to_write, [type_of_data, '_R.', resample_mesh, '.', fsLR_extension]));
lh_data_resample = lh_data_resample_gifti.cdata;
rh_data_resample = rh_data_resample_gifti.cdata;

%% Remove the output folder
rmdir(folder_to_write,'s');







