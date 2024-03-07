function  CBIG_MSHBM_example_single_subject(out_dir)

% CBIG_MSHBM_example_single_subject(out_dir)
%
% This function will call a simplified wrapper to run a single subject parcellation using MSHBM. This example use
% HCP trained group priors in fsaverage6 space.
% Results of the parcellation will be saved in the specified output directory. 
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
cd(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM'));

%% Prepare inputs
% prepare project_directory
params.project_dir = out_dir;

% prepare group prior
% user can set their own group prior
%params.group_prior = '<your_own>/priors/Params_Final.mat';

% prepare censor list
params.censor_list = {fullfile(CBIG_CODE_DIR,...
'/data/example_data/CoRR_HNU/subj01/subj01_sess1/qc',...
'subj01_sess1_bld002_FDRMS0.2_DVARS50_motion_outliers.txt')};
% prepare fMRI list
params.lh_fMRI_list = {fullfile(CBIG_CODE_DIR,...
'/data/example_data/CoRR_HNU/subj01/subj01_sess1/surf',...
'lh.subj01_sess1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz')};
params.rh_fMRI_list = {fullfile(CBIG_CODE_DIR,...
'/data/example_data/CoRR_HNU/subj01/subj01_sess1/surf',...
'rh.subj01_sess1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz')};

% data space
params.target_mesh = 'fsaverage6';
params.w = '50';
params.c = '50';
[lh_labels, rh_labels] = CBIG_MSHBM_parcellation_single_subject(params);