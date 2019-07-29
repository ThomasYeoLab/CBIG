function [lh_FS7_data,rh_FS7_data] = CBIG_project_fsLR2fsaverage(lh_fsLR_data,rh_fsLR_data,fsLR_mesh,type_of_data,folder_to_write,registration_version)

% [lh_FS7_data,rh_FS7_data] = CBIG_project_fsLR2fsaverage(lh_fsLR_data,rh_fsLR_data,fsLR_mesh,type_of_data,folder_to_write,registration_version)
%
% This function projects label/metric data in fs_LR_32k/fs_LR_164k to 
% fsaverage. The projection is performed in fs_LR_164k, therefore, if data 
% is in fs_LR_32k, the data will upsample to fs_LR_164k.
% Input:
%      -lh_fsLR_data, rh_fsLR_data:
%       a N x 1 and 1 x N vector, where N is the number of vertices. ?h_fsLR_data can
%       be metric data (float) or label data (integer).
%
%      -fsLR_mesh:
%       'fs_LR_32k'/'fs_LR_164k'. Mesh name of lh_fsLR_data, rh_fsLR_data.
%     
%      -type_of_data:
%       'metric' or 'label'. If ?h.fsLR_data is a metric data (float), then
%       type_of_data should be set to 'metric'. IF ?h.fsLR_data is a label 
%       data (integer), then type_of_data should be set to 'label'.
%
%      -folder_to_write:
%       output path. e.g. '/data/Mapping_FS_fsLR'.
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
%      -lh_FS7_data, rh_FS7_data:
%       output data in fsaverage after projection.
%
% Example:
% folder_to_write='/data/users/rkong/storage/ruby/data/HCP_relevant/Mapping_fsLR2FS';
% [lh_FS7_data,rh_FS7_data]=CBIG_project_fsLR2fsaverage(lh_label_fsLR_32k,rh_label_fsLR_32k,'fs_LR_32k','label',folder_to_write);
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% This function does not need vector check because the function itself
% contains checking statement.

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

%% input data should be Nx1 vector, if not, rotate it to be Nx1
if(size(lh_fsLR_data,2) ~= 1);
    lh_fsLR_data = lh_fsLR_data';
end
if(size(rh_fsLR_data,2) ~= 1);
    rh_fsLR_data = rh_fsLR_data';
end 


%% Save input ?h_fsLR_data as a cifti file
% set gifti file extension
if(strcmp(type_of_data, 'label'))
    fsLR_extension = 'label.gii'; % if type_of_data is label, extension of gifti file is .label.gii
elseif(strcmp(type_of_data, 'metric'))% assuming float data
    fsLR_extension = 'func.gii'; % if type_of_data is metric, extension of gifti file is .func.gii
else
    error('unknown type of data')
end

% get cifti func file structure
lh_fsLR_target = gifti(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', fsLR_mesh, 'example_func', 'Parcels_L.func.gii')); 
rh_fsLR_target = gifti(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', fsLR_mesh, 'example_func', 'Parcels_L.func.gii')); 

lh_fsLR_target.cdata = lh_fsLR_data;
rh_fsLR_target.cdata = rh_fsLR_data;

% save as a gifti file
save(lh_fsLR_target, fullfile(folder_to_write, [type_of_data, '_L.',fsLR_mesh, '.', fsLR_extension]), 'Base64Binary');
save(rh_fsLR_target, fullfile(folder_to_write, [type_of_data, '_R.',fsLR_mesh, '.', fsLR_extension]), 'Base64Binary');

%% In CBIG_project_fsLR2fsaverage.sh, it will
%  1) If fsLR_mesh is 'fs_LR_32k', it will upsample to fs_LR_164k
%  2) project from fs_LR_164k to fsaverage
system([fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'fslr_matlab', 'CBIG_project_fsLR2fsaverage.sh'), ' ', folder_to_write, ' ', type_of_data, ' ', fsLR_mesh, ' ', registration_version]);

%% Output fsaverage data
lh_FS7 = gifti(fullfile(folder_to_write, [type_of_data, '_L_gifti.gii']));
rh_FS7 = gifti(fullfile(folder_to_write, [type_of_data, '_R_gifti.gii']));

lh_FS7_data = lh_FS7.cdata;
rh_FS7_data = rh_FS7.cdata;

%% Remove the output folder
rmdir(folder_to_write,'s');
