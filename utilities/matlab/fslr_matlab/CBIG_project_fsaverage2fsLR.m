function [lh_fsLR_32k_data,rh_fsLR_32k_data,lh_fsLR_164k_data,rh_fsLR_164k_data] = CBIG_project_fsaverage2fsLR(lh_FS_data,rh_FS_data,FS_mesh,type_of_data,folder_to_write,registration_version)

% [lh_fsLR_32k_data,rh_fsLR_32k_data,lh_fsLR_164k_data,rh_fsLR_164k_data] = CBIG_project_fsaverage2fsLR(lh_FS_data,rh_FS_data,FS_mesh,type_of_data,folder_to_write,registration_version)
%
% This function projects label/metric data in
% fsaverage5/fsaverage6/fsaverage to fs_LR_32k/fs_LR_164k. The projection 
% is performed in fsaverage, therefore, if data is in fsaverage5/6, the 
% data will upsample to fsaverage.
% Input:
%      -lh_FS_data, rh_FS_data:
%       a 1 x N vector, where N is the number of vertices. ?h_FS_data can
%       be metric data (float) or label data (integer). 
%
%      -FS_mesh:
%       'fsaverage5'/'fsaverage6'/'fsaverage'. Mesh name of lh_FS_data, 
%       rh_FS_data.
%     
%      -type_of_data:
%       'metric' or 'label'. If ?h.FS_data is a metric data (float), then
%       type_of_data should be set to 'metric'. IF ?h.FS_data is a label 
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
%      -lh_label_fsLR_32k, rh_label_fsLR_32k:
%       output data in fs_LR_32k after projection.
%
%      -lh_label_fsLR_164k, lh_label_fsLR_164k:
%       output data in fs_LR_164k after projection.
% Example:
% load('/mnt/eql/yeo1/data/PublicParcellations/Schaefer/400_parcels/fsaverage6/mat_files/mat_files.mat')
% lh_FS_data=lh_label;rh_FS_data=rh_label;
% folder_to_write='/data/users/rkong/storage/ruby/data/HCP_relevant/Mapping_FS2fsLR'
% [lh_label_fsLR_32k,rh_label_fsLR_32k,lh_label_fsLR_164k,rh_label_fsLR_164k]=CBIG_project_fsaverage2fsLR(lh_FS_data,rh_FS_data,'fsaverage6','label',folder_to_write);

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

%% input data should be 1xN vector. Nx1 vector is not allowed in MARS_NNInterpolate_kdTree, MARS_linearInterpolate_kdTree
if(size(lh_FS_data,1) ~= 1);
    lh_FS_data = lh_FS_data';
end
if(size(rh_FS_data,1) ~= 1);
    rh_FS_data = rh_FS_data';
end 

%% upsample from fsaverage5/6 to fsaverage
if (~strcmp(FS_mesh, 'fsaverage'))% if its fsaverage5/6 data, we upsample from fsaverag5/6 to fsaverage
    lh_mesh7 = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
    rh_mesh7 = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex');
    lh_FS_mesh = CBIG_ReadNCAvgMesh('lh', FS_mesh, 'sphere', 'cortex');
    rh_FS_mesh = CBIG_ReadNCAvgMesh('rh', FS_mesh, 'sphere', 'cortex');
    
    if(strcmp(type_of_data, 'label')) % if type_of_data is label
        lh_data7 = MARS_NNInterpolate_kdTree(lh_mesh7.vertices, lh_FS_mesh, lh_FS_data);
        rh_data7 = MARS_NNInterpolate_kdTree(rh_mesh7.vertices, rh_FS_mesh, rh_FS_data); 
    elseif(strcmp(type_of_data, 'metric')) % if type_of_data is label
        lh_data7 = MARS_linearInterpolate_kdTree(lh_mesh7.vertices, lh_FS_mesh, lh_FS_data);
        rh_data7 = MARS_linearInterpolate_kdTree(rh_mesh7.vertices, rh_FS_mesh, rh_FS_data); 
    else
        error('unknown type of data')
    end   
elseif (strcmp(FS_mesh,'fsaverage'))% if its fsaverage data, we don't need to do upsampling
    % we already have fsaverage data
    lh_data7 = lh_FS_data;
    rh_data7 = rh_FS_data;
else
    error('The projection needs to be perform on fsaverage5/fsaverage6/fsaverage')
end

%% write out input data as .mgh file, so that this .mgh file can be later converted to gifti format
if(strcmp(type_of_data, 'label'))
    lh_mri.vol = int32(lh_data7);
    rh_mri.vol = int32(rh_data7);
    err = MRIwrite(lh_mri, fullfile(folder_to_write, 'label_L_fsaverage_borders.mgh'), 'int');
    err = MRIwrite(rh_mri, fullfile(folder_to_write, 'label_R_fsaverage_borders.mgh'), 'int');
    out_extension = 'label.gii'; % if type_of_data is label, extension of output file is .label.gii
elseif(strcmp(type_of_data, 'metric'))% assuming float data
    warning('assuming double data')
    lh_mri.vol = double(lh_data7);
    rh_mri.vol = double(rh_data7);
    err = MRIwrite(lh_mri, fullfile(folder_to_write, 'metric_L_fsaverage_borders.mgh'), 'double');
    err = MRIwrite(rh_mri, fullfile(folder_to_write, 'metric_R_fsaverage_borders.mgh'), 'double');
    out_extension = 'func.gii'; % if type_of_data is metric, extension of output file is .func.gii
else
    error('unknown type of data')
end

%% In CBIG_project_fsaverage2fsLR.sh, it will
%  1) convert .mgh file to gifti format
%  2) project to fs_LR_164k
%  3) downsample to fs_LR_32k
system([fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'fslr_matlab', 'CBIG_project_fsaverage2fsLR.sh'), ' ', folder_to_write, ' ', type_of_data, ' ', registration_version]);

%% Output fs_LR_32k/fs_LR_164k data
lh_fsLR_164k_gifti = gifti(fullfile(folder_to_write, [type_of_data, '_L.fs_LR_164k.', out_extension]));
rh_fsLR_164k_gifti = gifti(fullfile(folder_to_write, [type_of_data, '_R.fs_LR_164k.', out_extension]));
lh_fsLR_164k_data = lh_fsLR_164k_gifti.cdata;
rh_fsLR_164k_data = rh_fsLR_164k_gifti.cdata;

lh_fsLR_32k_gifti = gifti(fullfile(folder_to_write, [type_of_data, '_L.fs_LR_32k.', out_extension]));
rh_fsLR_32k_gifti = gifti(fullfile(folder_to_write, [type_of_data, '_R.fs_LR_32k.', out_extension]));
lh_fsLR_32k_data = lh_fsLR_32k_gifti.cdata;
rh_fsLR_32k_data = rh_fsLR_32k_gifti.cdata;

%% Remove the output folder
rmdir(folder_to_write,'s');
