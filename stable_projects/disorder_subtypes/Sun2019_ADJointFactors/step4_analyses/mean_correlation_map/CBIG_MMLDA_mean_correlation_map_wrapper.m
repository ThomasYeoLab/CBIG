function CBIG_MMLDA_mean_correlation_map_wrapper(out_dir)
% CBIG_MMLDA_mean_correlation_map_wrapper(out_dir)
%
% Wrapper function to get mean correlation map between mean top5 behavioral zscores and  atrophy zscores.
%
% Input:
%   - out_dir   : output directory
%
% Example:
%   CBIG_MMLDA_mean_correlation_map_wrapper('~/example')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

subinfo_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/step2_MMLDA/data'];
proj_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];

%%%
% get the z score of ADNI1 and ADNIGO2 AD 
%%%
% We don't use fieldnames, so ADNI1 is the same as ADNI2
[ADAS, MMSE, NEUROBAT] = CBIG_MMLDA_behavior_choose('ADNI2');
behavior_descript = [ADAS.DESCRIPTION' MMSE.DESCRIPTION' NEUROBAT.DESCRIPTION'];
behavior_sign = [ADAS.SIGN' MMSE.SIGN' NEUROBAT.SIGN'];
behavior_cate = [ADAS.CATE' MMSE.CATE' NEUROBAT.SIGN'];

behavior_fieldname_comb = [behavior_descript(1:4), {'ADAS:Naming+MMSE:Language'}, behavior_descript(6:8), ...
    {'ADAS:Spoken Language+Comprehension'}, behavior_descript([11 13:15]), {'MMSE:Immediate Recall+Delayed Recall'}, ...
    behavior_descript([17 20:21]), {'BNT: without cue+semantic cue'}, behavior_descript(27:end)];
behavior_sign_comb = [behavior_sign(1:4), {'PLUS'}, behavior_sign(6:8), {'PLUS'}, behavior_sign([11 13:15]), ...
    {'MINUS'}, behavior_sign([17 20:21]), {'MINUS'}, behavior_sign(27:end)];
behavior_cate_comb = [behavior_cate(1:4), {'NONE'}, behavior_cate(6:8), {'NONE'}, behavior_cate([11 13:15]), ...
    {'MEM'}, behavior_cate([17 20:21]), {'NONE'}, behavior_cate(27:end)];

% get ADNIGO2 AD zscore
subinfo = csvread([subinfo_dir '/ADNI2_bl_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

% combine behaviors
behavior_comb = [behavior(:, 1:4), behavior(:, 5)+2-behavior(:, 19), behavior(:, 6:8), ...
    behavior(:, 10)+behavior(:, 12), behavior(:, [11 13:15]), behavior(:, 16)+behavior(:, 18), ...
    behavior(:, [17 20:21]), sum(behavior(:, 23:24), 2), behavior(:, 27:end)];

% remove behavior with missing scores > 5
num_of_nan = sum(isnan(behavior_comb), 2);
ind = (num_of_nan <= 5);
character_rm = character(ind, :);
behavior_comb_rm = behavior_comb(ind, :);

% matrix completion for missing scores
behavior_comb_rm_mc = behavior_comb_rm;
ind_AD = character_rm(:, 4) == 3;
ind_ALL_noNAN = ~isnan(sum(behavior_comb_rm, 2));
behavior_comb_rm_AD_mc = CBIG_MMLDA_matrix_completion_GLM(behavior_comb_rm(ind_AD, :), ...
    behavior_comb_rm(ind_ALL_noNAN, :));
behavior_comb_rm_mc(ind_AD, :) = behavior_comb_rm_AD_mc;

% remove behavior with missing scores
ind_noNAN = ~isnan(sum(behavior_comb_rm_mc, 2));
character_rm_rm = character_rm(ind_noNAN, :);
behavior_comb_rm_mc_rm = behavior_comb_rm_mc(ind_noNAN, :);

% Brain: brain to zscore
ind_brain = CBIG_MMLDA_find_array_in_array(character_rm_rm(:, 1), character(:, 1));
vol4D = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI2_bl/mri/gm_merg_s10_MNI2mm.nii.gz']);
vol4D.vol = vol4D.vol(:, :, :, ind_brain);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz']);
regressor = character(ind_brain, [3 2 5]);
refpara = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/ADNIDataRelease/regression_zscore_paras/ADNI2_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb_rm_mc_rm;
regressor = character_rm_rm(:, [3 2 5]);
refpara = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/ADNIDataRelease/regression_zscore_paras/ADNI2_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

ADNIGO2_AD_Z_brain = Z_brain(character_rm_rm(:, 4) == 3, :);
ADNIGO2_AD_Z_behavior = Z_behavior(character_rm_rm(:, 4) == 3, :);

% get ADNI1 AD zscore
subinfo = csvread([subinfo_dir '/ADNI1_bl_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

% combine behaviors
behavior_comb = [behavior(:, 1:4), behavior(:, 5)+2-behavior(:, 19), behavior(:, 6:8), ...
    behavior(:, 10)+behavior(:, 12), behavior(:, [11 13:15]), behavior(:, 16)+behavior(:, 18), ...
    behavior(:, [17 20:21]), sum(behavior(:, 23:24), 2), behavior(:, 27:end)];

% remove behavior with missing scores
ind_noNAN = ~isnan(sum(behavior_comb, 2));
character_rm = character(ind_noNAN, :);
behavior_comb_rm = behavior_comb(ind_noNAN, :);

% Brain: brain to zscore
ind_brain = CBIG_MMLDA_find_array_in_array(character_rm(:, 1), character(:, 1));
vol4D = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI1_bl/mri/gm_merg_s10_MNI2mm.nii.gz']);
vol4D.vol = vol4D.vol(:, :, :, ind_brain);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNI1/ADNI1_bl_gm_mask_MNI2mm.nii.gz']);
regressor = character(ind_brain, [3 2 5]);
refpara = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/ADNIDataRelease/regression_zscore_paras/ADNI1_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb_rm;
regressor = character_rm(:, [3 2 5]);
refpara = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/ADNIDataRelease/regression_zscore_paras/ADNI1_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

ADNI1_AD_Z_brain = Z_brain(character_rm(:, 4) == 3, :);
ADNI1_AD_Z_behavior = Z_behavior(character_rm(:, 4) == 3, :);

% get the index of top5 scores
beta2 = exp(load([proj_dir '/step2_MMLDA/results/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1/k3/r9/final.beta2']));
beta2_order = beta2([1 3 2], :);
[~, ind_MEM] = sort(beta2_order(1, :), 'descend');
[~, ind_LAN] = sort(beta2_order(2, :), 'descend');
[~, ind_EF]  = sort(beta2_order(3, :), 'descend');
MEM_top5 = ind_MEM(1:5);
LAN_top5 = ind_LAN(1:5);
EF_top5  = ind_EF(1:5);

% compute mean of top5 zscores and correlate it with atrophy maps
ADNIGO2_MEM = CBIG_corr(mean(ADNIGO2_AD_Z_behavior(:, MEM_top5), 2), ADNIGO2_AD_Z_brain);
ADNIGO2_LAN = CBIG_corr(mean(ADNIGO2_AD_Z_behavior(:, LAN_top5), 2), ADNIGO2_AD_Z_brain);
ADNIGO2_EF = CBIG_corr(mean(ADNIGO2_AD_Z_behavior(:, EF_top5), 2), ADNIGO2_AD_Z_brain);
ADNI1_MEM = CBIG_corr(mean(ADNI1_AD_Z_behavior(:, MEM_top5), 2), ADNI1_AD_Z_brain);
ADNI1_LAN = CBIG_corr(mean(ADNI1_AD_Z_behavior(:, LAN_top5), 2), ADNI1_AD_Z_brain);
ADNI1_EF = CBIG_corr(mean(ADNI1_AD_Z_behavior(:, EF_top5), 2), ADNI1_AD_Z_brain);

%%%
% visualize mean correlation map in MNI space
%%%
underlay_vol = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA' ...
    '/Template_T1_IXI555_MNI152_brain_MNI2mm.nii'];
color_map_name = 'ColdHot';
plane = 'coronal';
out_dir = [out_dir '/mean_correlation_map'];
mkdir(out_dir)

% ADNIGO2
mask = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/ADNIDataRelease' ...
    '/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz'];
in_vol = [out_dir '/ADNIGO2_MEM.nii.gz'];
CBIG_MMLDA_vector2nifti(ADNIGO2_MEM, mask, in_vol);
min_thresh = 0.1; 
max_thresh = 0.3;
out_name = 'ADNIGO2_AD_MEM';
CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
    color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)

in_vol = [out_dir '/ADNIGO2_LAN.nii.gz'];
CBIG_MMLDA_vector2nifti(ADNIGO2_LAN, mask, in_vol);
min_thresh = 0.2; 
max_thresh = 0.4;
out_name = 'ADNIGO2_AD_LAN';
CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
    color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)

in_vol = [out_dir '/ADNIGO2_EF.nii.gz'];
CBIG_MMLDA_vector2nifti(ADNIGO2_EF, mask, in_vol);
min_thresh = 0.2; 
max_thresh = 0.4;
out_name = 'ADNIGO2_AD_EF';
CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
    color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)

% ADNI1
mask = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/ADNIDataRelease' ...
    '/SPM_VBM_files/ADNI1/ADNI1_bl_gm_mask_MNI2mm.nii.gz'];
in_vol = [out_dir '/ADNI1_MEM.nii.gz'];
CBIG_MMLDA_vector2nifti(ADNI1_MEM, mask, in_vol);
min_thresh = 0.2; 
max_thresh = 0.4;
out_name = 'ADNI1_AD_MEM';
CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
    color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)

in_vol = [out_dir '/ADNI1_LAN.nii.gz'];
CBIG_MMLDA_vector2nifti(ADNI1_LAN, mask, in_vol);
min_thresh = 0.2; 
max_thresh = 0.4;
out_name = 'ADNI1_AD_LAN';
CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
    color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)

in_vol = [out_dir '/ADNI1_EF.nii.gz'];
CBIG_MMLDA_vector2nifti(ADNI1_EF, mask, in_vol);
min_thresh = 0.2; 
max_thresh = 0.4;
out_name = 'ADNI1_AD_EF';
CBIG_MMLDA_view_vol_slice_with_underlay(in_vol, underlay_vol, ...
    color_map_name, plane, min_thresh, max_thresh, out_dir, out_name)


rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])