function [lh_labels, rh_labels] = CBIG_ArealMSHBM_parcellation_single_subject(params)
% [lh_labels, rh_labels] = CBIG_ArealMSHBM_parcellation_single_subject(params)
%
% This function is a simplied wrapper to perform MSHBM individual areal-level
% parcellation for a single subject. This function prepare take simple inputs
% and prepare the input argumemts and required data lists for 
% CBIG_MSHBM_generate_profiles.m and CBIG_d/c/gMSHBM_generate_individual_parcellation.m.
% Since only gMSHBM group prior is available, this function only supports gMSHBM
% model if user wants to use default group priors. If user has their own group
% prior, they can specify the full path to the group prior file in the group_prior
% field in the input params and use cMSHBM or dMSHBM model.
%
% Inputs:
%   - params:
%     A structure with the following fields:
%       - project_dir:
%         A string. The directory to store the parcellation results. The project_dir
%         should be unique for each subject.
%
%       - censor_list:
%         A 1XN or Nx1 cell array of strings. Each string is the full path to the
%         censor file for each session. If there is only one session, censor_list
%         can be a string. For censor file format, it should be a txt file with
%         one column. Each row corresponds to a time point. 1 means the time point
%         is not censored, 0 means the time point is censored. If there is no
%         censor file or you don't want to censor any time point, you can set censor
%         list to be a file with 1 in each row.
%
%       - lh_fMRI_list:
%         A 1XN or Nx1 cell array of strings. Each string is the full path to the
%         fMRI file for each session. If there is only one session, lh_fMRI_list
%         can be a string.
%
%       - rh_fMRI_list: (should be empty or not specified if target_mesh is fs_LR_32k)
%         A 1XN or Nx1 cell array of strings. Each string is the full path to the
%         fMRI file for each session. If there is only one session, rh_fMRI_list
%         can be a string. If target_mesh is fs_LR_32k, rh_fMRI_list can be empty.
%
%       - target_mesh:
%         A string. The target mesh. It can be 'fsaverage6', or 'fs_LR_32k'. The
%         target_mesh can also be 'fsaverage5', but the default group prior is only
%         available in fsaverage6 and fs_LR_32k. If you want to use fsaverage5 as
%         the target mesh, you need to generate your own group prior and specify
%         the full path to the group prior file in the group_prior field.
%
%       - num_ROIs:
%         A string. The number of ROIs. Default is '400'.
%
%       - model:
%         A string. The model to use. 'gMSHBM, 'cMSHBM' or 'dMSHBM'. Default is 'gMSHBM'.
%         Currently only gMSHBM model is supported if user wants to use default group
%         priors. If user has their own group prior, they can specify the full path
%         to the group prior file in the group_prior field and use cMSHBM or dMSHBM model.
%
%       - group_prior:
%         A string. The full path to the group prior file "Params_Final.mat". If
%         group_prior is not specified, i.e., group_prior field doesn't exist in
%         params, the default group prior will be used. If target_mesh is 
%         fsaverage6, the default group prior is:
%         $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/
%         group_priors/HCP_fsaverage6_40sub/<num_ROIs>/gMSHBM/beta5/Params_Final.mat
%         If target_mesh is fs_LR_32k, the default group prior is:
%         $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/lib/
%         group_priors/HCP_fs_LR_32k_40sub/<num_ROIs>/gMSHBM/beta50/Params_Final.mat
%         If you generate your own group prior, you can specify the full path to
%         your group prior file to use it. Please note that we do not provide group
%         priors for fsaverage5. If you want to use fsaverage5 as the target mesh,
%         you need to generate your own group prior and specify the full path to the
%         group prior file in the group_prior field.
%
%       - beta: (no need to specify if model is dMSHBM)
%         A string. The weight of spatial localization prior. Default is '5' for
%         fsaverage6 and '50' for fs_LR_32k. Stronger weight will result in 
%         parcellation with more spatially localized ROIs.
%
%       - w:
%         A string. The weight for the group spatial prior. Default is '50'.
%         Stronger weight will result in parcellation spatially closer to the
%         group-level networks.
%
%       - c:
%         A string. The weight for the MRF constraint. Default is '10'. Stronger
%         weight will result in smoother parcellation.
%
%       - overwrite_flag:
%         A number. If overwrite_flag is 1, the project_dir will be deleted if it
%         already exists. Default is 0. Please note that if profiles already exist
%         in the project_dir, setting overwrite_flag to 1 will delete the existing
%         profiles and re-generate them.
%
% Outputs:
%   Results will be saved in project_dir:
%   <project_dir>/ind_parcellation_<model>/test_set
%   - lh_labels, rh_labels:
%     parcellation results for left and right hemisphere.
%
% Example for a subject with 1 run of fsaverage6 fMRI data:
%   params.project_dir = '/myproject/sub1';
%   params.censor_list = '/mydata/sub1/censor.txt';
%   params.lh_fMRI_list = '/mydata/sub1/lh.fsaverage6_surf.nii.gz';
%   params.rh_fMRI_list = '/mydata/sub1/rh.fsaverage6_surf.nii.gz';
%   params.target_mesh = 'fsaverage6';
%   CBIG_ArealMSHBM_parcellation_single_subject(params);
%
% Example for a subject with 2 runs of fs_LR_32k fMRI data:
%   params.project_dir = '/myproject/sub1';
%   params.censor_list = {'/mydata/sub1/censor1.txt', '/mydata/sub1/censor2.txt'};
%   params.lh_fMRI_list = {'/mydata/sub1/fs_LR_32k_sess1_surf.dtseries.nii',...
%  '/mydata/sub1/fs_LR_32k_sess2_surf.dtseries.nii'};
%   params.target_mesh = 'fs_LR_32k';
%   CBIG_ArealMSHBM_parcellation_single_subject(params);
%
% Example for a subject with 2 runs of fsaverage6 fMRI data and specify all parameters:
%   params.project_dir = '/myproject/sub1';
%   params.censor_list = {'/mydata/sub1/censor1.txt', '/mydata/sub1/censor2.txt'};
%   params.lh_fMRI_list = {'/mydata/sub1/lh.fsaverage6_sess1_surf.nii.gz',...
%  '/mydata/sub1/lh.fsaverage6_sess2_surf.nii.gz'};
%   params.rh_fMRI_list = {'/mydata/sub1/rh.fsaverage6_sess1_surf.nii.gz',...
%  '/mydata/sub1/rh.fsaverage6_sess2_surf.nii.gz'};
%   params.target_mesh = 'fsaverage6';
%   params.model = 'dMSHBM';
%   params.num_ROIs = '300';
%   params.w = '50';
%   params.c = '30';
%   params.group_prior = '/myproject/training_step/priors/dMSHBM/Params_Final.mat';
%   params.overwrite_flag = 1;
%   CBIG_ArealMSHBM_parcellation_single_subject(params);
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation','Kong2019_MSHBM',...
    'step1_generate_profiles_and_ini_params'));
addpath(genpath(fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation','Kong2022_ArealMSHBM',...
    'step0_generate_gradient_prior')));
addpath(genpath(fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation','Kong2022_ArealMSHBM',...
    'step3_generate_ind_parcellations')));

%% Prepare inputs
% prepare project_directory
if ~isfield(params, 'project_dir')
    error('project_dir is not specified!');
else
    project_dir = params.project_dir;
end
% prepare group prior
if isfield(params, 'group_prior')
    group_prior = params.group_prior;
end
% prepare target mesh
if ~isfield(params, 'target_mesh')
    target_mesh = 'fsaverage6';
else
    target_mesh = params.target_mesh;
end
% prepare censor list
if ~isfield(params, 'censor_list')
    error('censor_list is not specified!');
else
    censor_list = params.censor_list;
end
% prepare fMRI list
if ~isfield(params, 'lh_fMRI_list')
    error('lh_fMRI_list is not specified!');
else
    lh_fMRI_list = params.lh_fMRI_list;
end
if ~isfield(params, 'rh_fMRI_list')
    if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
        error('rh_fMRI_list is not specified!');
    else
        rh_fMRI_list = [];
    end
else
    rh_fMRI_list = params.rh_fMRI_list;
end
if ~isfield(params, 'model')
    model = 'gMSHBM';
else
    model = params.model;
end
if ~isfield(params, 'num_ROIs')
    num_ROIs = '400';
else
    num_ROIs = params.num_ROIs;
end
if ~isfield(params, 'beta')
    if strcmp(target_mesh, 'fsaverage6')
        beta = '5';
    elseif strcmp(target_mesh, 'fs_LR_32k')
        beta = '50';
    end
else
    beta = params.beta;
end
if ~isfield(params, 'w')
    w = '50';
else
    w = params.w;
end
if ~isfield(params, 'c')
    c = '10';
else
    c = params.c;
end
if ~isfield(params, 'overwrite_flag')
    overwrite_flag = 0;
else
    overwrite_flag = params.overwrite_flag;
end

%% Check inputs
% if overwrite_flag is not specified, use 0 as default
if ~exist('overwrite_flag', 'var')
    overwrite_flag = 0;
end
% if target_mesh is not specified, use fsaverage6
if ~exist('target_mesh', 'var')
    target_mesh = 'fsaverage6';
end
% if target mesh is fsaverage6 or fsaverage5, seed_mesh should be fsaverage3
if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
    seed_mesh = 'fsaverage3';
    % censor list and fmri list cannot be empty
    if isempty(censor_list) || isempty(lh_fMRI_list) || isempty(rh_fMRI_list)
        error('censor list and fMRI list cannot be empty!');
    end
% if target mesh is fs_LR_32k, seed_mesh should be fs_LR_900
elseif strcmp(target_mesh, 'fs_LR_32k')
    seed_mesh = 'fs_LR_900';
    % censor list cannot be empty
    if isempty(censor_list) || isempty(lh_fMRI_list)
        error('censor list and fMRI list cannot be empty!');
    end
end

% overwrite_flag is 1, delete the project_dir
if overwrite_flag == 1
    system(['rm -rf ' project_dir]);
end

% if group_prior is not specified, use the default one
% if target_mesh is fsaverage6, use:
% $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM...
%       /lib/group_priors/HCP_fsaverage6_40sub/<num_ROIs>/gMSHBM/beta<?>/Params_Final.mat
% if target_mesh is fs_LR_32k, use:
% $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM...
%       /lib/group_priors/HCP_fs_LR_32K_40sub/<num_ROIs>/gMSHBM/beta<?>/Params_Final.mat
% fsaverage5 group prior is not provided. If you want to use fsaverage5 as the
% target mesh, you need to generate your own group prior and specify the full
% path to the group prior file in the group_prior field.
beta_folder = ['beta' beta];
if ~exist('group_prior', 'var')
    if strcmp(target_mesh, 'fsaverage6')
        group_prior = fullfile(CBIG_CODE_DIR, 'stable_projects','brain_parcellation',...
        'Kong2022_ArealMSHBM','lib','group_priors','HCP_fsaverage6_40sub',...
         num_ROIs, 'gMSHBM', beta_folder,'Params_Final.mat');
    elseif strcmp(target_mesh, 'fs_LR_32k')
        group_prior = fullfile(CBIG_CODE_DIR, 'stable_projects','brain_parcellation',...
        'Kong2022_ArealMSHBM','lib','group_priors','HCP_fs_LR_32k_40sub',...
         num_ROIs, 'gMSHBM', beta_folder,'Params_Final.mat');
    else
        error('Only fsaverage6 and fs_LR_32k group priors are provided. Specify group_prior for fsaverage5!');
    end
end
% create a symbolic link to group_prior
if strcmp(model, 'cMSHBM') || strcmp(model, 'gMSHBM')
    if ~exist(fullfile(project_dir, 'priors', model, beta_folder), 'dir')
        mkdir(fullfile(project_dir, 'priors',model, beta_folder));
    end
    if ~exist(fullfile(project_dir, 'priors',model, beta_folder ,'Params_Final.mat'), 'file')
        system(['ln -s ' group_prior ' ' fullfile(project_dir, 'priors',model, beta_folder, 'Params_Final.mat')]);
    else
        disp('group_prior already exists!');
    end
elseif strcmp(model, 'dMSHBM')
    if ~exist(fullfile(project_dir, 'priors', model), 'dir')
        mkdir(fullfile(project_dir, 'priors',model));
    end
    if ~exist(fullfile(project_dir, 'priors',model ,'Params_Final.mat'), 'file')
        system(['ln -s ' group_prior ' ' fullfile(project_dir, 'priors',model, 'Params_Final.mat')]);
    else
        disp('group_prior already exists!');
    end
else
    error('model should be gMSHBM, cMSHBM or dMSHBM!');
end

% create symbolic link to spatial mask
if(exist(fullfile(project_dir,'spatial_mask'),'dir'))
    system(['rm -r ' fullfile(project_dir,'spatial_mask')]);
end
mkdir(fullfile(project_dir, 'spatial_mask'));
if strcmp(target_mesh, 'fsaverage6')
    spatial_mask = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation',...
        'Kong2022_ArealMSHBM','lib','spatial_mask', num_ROIs, 'spatial_mask_fsaverage6.mat');
elseif strcmp(target_mesh, 'fs_LR_32k')
    spatial_mask = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation',...
        'Kong2022_ArealMSHBM','lib','spatial_mask', num_ROIs, 'spatial_mask_fs_LR_32k.mat');
end
system(['ln -s ' spatial_mask ' ' fullfile(project_dir, 'spatial_mask', ['spatial_mask_' target_mesh '.mat'])]);

% censor_list and fMRI_list should be a 1XN cell array or Nx1 cell array
% check if censor_list is a cell array
censor_flag = 1;
if iscell(censor_list)
    if min(size(censor_list)) ~= 1
        error('censor_list should be a 1XN or Nx1 cell array!');
    end
    % check number of sessions
    num_sess = length(censor_list);
elseif ischar(censor_list)
    if strcmp(censor_list, 'NONE')
        censor_flag = 0;
        num_sess = length(lh_fMRI_list);
    else
        censor_flag = 1;
        num_sess = 1;
    end
else
    error('censor_list should be a 1XN or Nx1 cell array!');
end
% check if fMRI_list is a cell array
if iscell(lh_fMRI_list)
    if min(size(lh_fMRI_list)) ~= 1
        error('lh_fMRI_list should be a 1XN or Nx1 cell array!');
    end
elseif ischar(lh_fMRI_list)
    lh_fMRI_list = {lh_fMRI_list};
end
if iscell(rh_fMRI_list)
    if min(size(rh_fMRI_list)) ~= 1
        error('rh_fMRI_list should be a 1XN or Nx1 cell array!');
    end
elseif ischar(rh_fMRI_list)
    rh_fMRI_list = {rh_fMRI_list};
end

% if target mesh is fsaverage6 or fsaverage5, 
% check if number of sessions are equal between censor list and fMRI list
if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
    if num_sess ~= length(lh_fMRI_list) || num_sess ~= length(rh_fMRI_list)
        error('Number of sessions are not equal between censor list and fMRI list!');
    end
end
% if target mesh is fs_LR_32k, 
% check if number of sessions are equal between censor list and lh_fMRI_list or rh_fMRI_list
if strcmp(target_mesh, 'fs_LR_32k')
    if num_sess ~= length(lh_fMRI_list) && num_sess ~= length(rh_fMRI_list)
        error('Number of sessions are not equal between censor list and fMRI list!');
    end
end    

% export data list
if censor_flag == 1
    censor_list_filedir = fullfile(project_dir, 'data_list', 'censor_list');
    save_cell_as_list(censor_list, censor_list_filedir, 'sub1_sess');
end
fMRI_list_filedir = fullfile(project_dir, 'data_list', 'fMRI_list');

% if target mesh is fsaverage6 or fsaverage5, export lh_fMRI_list and rh_fMRI_list
if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
    save_cell_as_list(lh_fMRI_list, fMRI_list_filedir, 'lh_sub1_sess');
    save_cell_as_list(rh_fMRI_list, fMRI_list_filedir, 'rh_sub1_sess');
end
% if target mesh is fs_LR_32k, export the lh_fMRI_list
if strcmp(target_mesh, 'fs_LR_32k')
    save_cell_as_list(lh_fMRI_list, fMRI_list_filedir, 'sub1_sess');
end

if strcmp(model, 'gMSHBM')
        
    % prepare gradients
    do_gradient = 1;
    gradient_save_dir = fullfile(project_dir, 'gradients','sub1');
    gradient_name = '_emb_100_distance_matrix.mat';
    lh_gradient_list = fullfile(gradient_save_dir, ['lh' gradient_name]);
    rh_gradient_list = fullfile(gradient_save_dir, ['rh' gradient_name]);
    if exist(lh_gradient_list, 'file') && exist(rh_gradient_list, 'file')
        do_gradient = 0;
        disp('lh and rh gradients already exist!');
    end
    if do_gradient == 1
        CBIG_ArealMSHBM_generate_gradient(target_mesh, project_dir, '1' , num2str(num_sess));
    end
    % export gradient list
    gradient_list_filedir = fullfile(project_dir, 'gradient_list','test_set');
    save_cell_as_list(lh_gradient_list, gradient_list_filedir, 'gradient_list_lh', 0);
    save_cell_as_list(rh_gradient_list, gradient_list_filedir, 'gradient_list_rh', 0);
end

% prepare profiles
do_profile = 1;
if num_sess == 1
    profile_save_dir = fullfile(project_dir, 'profiles', 'sub1', 'sess1');
    profile_name_fake1 = ['sub1_sess1_' target_mesh '_roi' seed_mesh '.surf2surf_profile_1'];
    profile_name_fake2 = ['sub1_sess1_' target_mesh '_roi' seed_mesh '.surf2surf_profile_2'];
    % for target mesh fsaverage6 or fsaverage5, save lh and rh profiles
    if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
        lh_profile_list = {fullfile(profile_save_dir, ['lh.' profile_name_fake1 '.nii.gz']),...
            fullfile(profile_save_dir, ['lh.' profile_name_fake2 '.nii.gz'])};
        rh_profile_list = {fullfile(profile_save_dir, ['rh.' profile_name_fake1 '.nii.gz']),...
            fullfile(profile_save_dir, ['rh.' profile_name_fake2 '.nii.gz'])};
        % check if lh and rh profiles exist
        if exist(lh_profile_list{1}, 'file') && exist(rh_profile_list{1}, 'file')
            do_profile = 0;
            disp('lh and rh profiles already exist!');
        end
    end
    % for target mesh fs_LR_32k, save profiles
    if strcmp(target_mesh, 'fs_LR_32k')
        profile_list = {fullfile(profile_save_dir, [profile_name_fake1 '.mat']),...
            fullfile(profile_save_dir, [profile_name_fake2 '.mat'])};
        % check if profiles exist
        if exist(profile_list{1}, 'file')
            do_profile = 0;
            disp('profiles already exist!');
        end
    end
    if do_profile == 1
        CBIG_MSHBM_generate_profiles(seed_mesh, target_mesh, project_dir, '1', '1', '1');
    end

else
    for curr_sess = 1:num_sess
        profile_save_dir = fullfile(project_dir, 'profiles', 'sub1', ['sess' num2str(curr_sess)]);
        profile_name = ['sub1_sess' num2str(curr_sess) '_' target_mesh '_roi' seed_mesh '.surf2surf_profile'];
        % for target mesh fsaverage6 or fsaverage5, save lh and rh profiles
        if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
            lh_profile_list{curr_sess} = fullfile(profile_save_dir, ['lh.' profile_name '.nii.gz']);
            rh_profile_list{curr_sess} = fullfile(profile_save_dir, ['rh.' profile_name '.nii.gz']);
            % check if lh and rh profiles exist
            if exist(lh_profile_list{curr_sess}, 'file') && exist(rh_profile_list{curr_sess}, 'file')
                do_profile = 0;
                disp('lh and rh profiles already exist!');
            end
        end
        % for target mesh fs_LR_32k, save profiles
        if strcmp(target_mesh, 'fs_LR_32k')
            profile_list{curr_sess} = fullfile(profile_save_dir, [profile_name '.mat']);
            % check if profiles exist
            if exist(profile_list{curr_sess}, 'file')
                do_profile = 0;
                disp('profiles already exist!');
            end
        end
        if do_profile == 1
            CBIG_MSHBM_generate_profiles(seed_mesh, target_mesh, project_dir, '1', num2str(curr_sess), '0');
        end
    end
end

% export profile list
profile_list_filedir = fullfile(project_dir, 'profile_list','test_set');
% for target mesh fsaverage6 or fsaverage5, export lh and rh profiles
if strcmp(target_mesh, 'fsaverage6') || strcmp(target_mesh, 'fsaverage5')
    save_cell_as_list(lh_profile_list, profile_list_filedir, 'lh_sess');
    save_cell_as_list(rh_profile_list, profile_list_filedir, 'rh_sess');
end
% for target mesh fs_LR_32k, export profiles
if strcmp(target_mesh, 'fs_LR_32k')
    save_cell_as_list(profile_list, profile_list_filedir, 'sess');
end

% perform parcellation
if(num_sess == 1)
    if strcmp(model, 'gMSHBM')
        [lh_labels, rh_labels] = CBIG_ArealMSHBM_gMSHBM_generate_individual_parcellation(project_dir,...
            target_mesh, '2', num_ROIs, '1', w, c, beta, 'test_set');
    elseif strcmp(model, 'cMSHBM')
        [lh_labels, rh_labels] = CBIG_ArealMSHBM_cMSHBM_generate_individual_parcellation(project_dir,...
            target_mesh, '2', num_ROIs, '1', w, c, beta, 'test_set');
    elseif strcmp(model, 'dMSHBM')
        [lh_labels, rh_labels] = CBIG_ArealMSHBM_dMSHBM_generate_individual_parcellation(project_dir,...
            target_mesh, '2', num_ROIs, '1', w, c, 'test_set');
    end
else
    [lh_labels, rh_labels] = CBIG_ArealMSHBM_gMSHBM_generate_individual_parcellation(project_dir,...
        target_mesh, num2str(num_sess), num_ROIs, '1', w, c, beta, 'test_set');
end

rmpath(fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation','Kong2019_MSHBM',...
    'step1_generate_profiles_and_ini_params'));
rmpath(genpath(fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation','Kong2022_ArealMSHBM',...
    'step0_generate_gradient_prior')));
rmpath(genpath(fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation','Kong2022_ArealMSHBM',...
    'step3_generate_ind_parcellations')));

end

function save_cell_as_list(myCell, filedir, baseFileName, idx_flag)
% save_cell_as_list(myCell, filedir, baseFileName, idx_flag)
%
% This function saves a cell array as a list of files.
% Each row in the cell array will be saved in a separate file.
% if myCell is a string, convert it to cell
% For a single row cell, if idx_flag is 0, do not add index to the file name
if (nargin < 4)
    idx_flag = 1;
end
if ischar(myCell)
    myCell = {myCell};
end

% check if filedir exists
if ~exist(filedir, 'dir')
    mkdir(filedir);
end

if (length(myCell) == 1 && idx_flag == 0)
    % Open the file for writing
    fileID = fopen(fullfile(filedir, [baseFileName '.txt']), 'w');
    
    % Check if the file was opened successfully
    if fileID == -1
        error('Unable to open file for writing: %s', fullfile(filedir, [baseFileName '.txt']));
    end
    
    % Write the cell content to the file
    fprintf(fileID, '%s\n', myCell{1});
    
    % Close the file
    fclose(fileID);
    
    disp(['Row 1 has been successfully exported to ' fullfile(filedir, [baseFileName '.txt'])]);
else
    % Iterate through the cell array and save each row in a separate file
    for row = 1:length(myCell)
        % Generate the filename with the row index
        fileName = sprintf(fullfile(filedir, '%s%d.txt'), baseFileName, row);
        
        % Open the file for writing
        fileID = fopen(fileName, 'w');
        
        % Check if the file was opened successfully
        if fileID == -1
            error('Unable to open file for writing: %s', fileName);
        end
        
        % Write the cell content to the file
        fprintf(fileID, '%s\n', myCell{row});
        
        % Close the file
        fclose(fileID);
        
        disp(['Row ' num2str(row) ' has been successfully exported to ' fileName]);
    end
end
    

end