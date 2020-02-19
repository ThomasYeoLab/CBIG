function CBIG_MMLDA_brain_behavior_to_doc_wrapper(subinfo_dir, out_dir)
% This wrapper converts brain nifti files and behavior to documents used for MMLDA estimation.
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 1 || isempty(subinfo_dir)
    subinfo_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/step2_MMLDA/data'];
end
if nargin < 2 || isempty(out_dir)
    out_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/step2_MMLDA/results/BrainBehavior2doc'];
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ADNI_VBM_path = [getenv('CBIG_MMLDA_ANDI_DATA') '/Sun2019_SPMVBM'];

addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

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


%%%
% ADNI2 bl 
%%%
subinfo = csvread([subinfo_dir '/ADNI2_bl_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

% remove ADAS: Recall instruction because CN subjects has 0 std
% remove ANART because it doesn't exist in m6 m12 ...
% combine BNT:without cue and semantic cue
% remove CFT: vegetable because it doesn't exist in ADNI2
% combine MMSE: Language and ADAS: Naming
% combine MMSE: Immediate Recall with MMSE: Delayed Recall
% combine ADAS: Spoken Language, Comprehension
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

% Brain: create reference paramters 
vol4D_bl = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI2_bl/mri/gm_merg_MNI2mm_s10.nii.gz']);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz']);
ind_CN_bl_brain = (subinfo(:, 4) == 1);
regressor_bl_brain = subinfo(:, [3 2 5]);
zscore_method = 'meanCNstdALL';
output = [out_dir '/ADNI2_bl_brain_refpara.mat'];
CBIG_MMLDA_brain_to_zscore_create_refpara(vol4D_bl, mask, ind_CN_bl_brain, regressor_bl_brain, zscore_method, output);

% Behavior: create reference paramters 
behavior_bl = behavior_comb_rm_mc_rm;
ind_CN_bl_behavior = (character_rm_rm(:, 4) == 1);
regressor_bl_behavior = character_rm_rm(:, [3 2 5]);
zscore_method = 'meanCNstdALL';
output = [out_dir '/ADNI2_bl_behavior_refpara.mat'];
CBIG_MMLDA_behavior_to_zscore_create_refpara(behavior_bl, ind_CN_bl_behavior, ...
    regressor_bl_behavior, zscore_method, output);

% Brain: brain to zscore
ind_brain = CBIG_MMLDA_find_array_in_array(character_rm_rm(:, 1), character(:, 1));
vol4D = vol4D_bl;
vol4D.vol = vol4D.vol(:, :, :, ind_brain);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz']);
regressor = character(ind_brain, [3 2 5]);
refpara = [out_dir '/ADNI2_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb_rm_mc_rm;
regressor = character_rm_rm(:, [3 2 5]);
refpara = [out_dir '/ADNI2_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

% write zscore to doc
rid = character_rm_rm(:, 1);
ind_out = (character_rm_rm(:, 4) == 3);
plus_num = 1;
nfold = 10;
out_name = 'ADNI2_bl_AD_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

rid = character_rm_rm(:, 1);
ind_out = (character_rm_rm(:, 4) == 3);
plus_num = 1;
nfold = 1;
out_name = 'ADNI2_bl_AD_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

rid = character_rm_rm(:, 1);
ind_out = (character_rm_rm(:, 4) == 2);
plus_num = 1;
nfold = 1;
out_name = 'ADNI2_bl_MCI_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

rid = character_rm_rm(:, 1);
ind_out = (character_rm_rm(:, 4) == 1);
plus_num = 1;
nfold = 1;
out_name = 'ADNI2_bl_CN_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

%%%
% ADNI2 m12
%%%
% Load ADNI2 bl subinfo to make sure m12 rid is from bl rid
subinfo = csvread([subinfo_dir '/ADNI2_bl_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

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
character_rm_rm_ADNI2_bl = character_rm(ind_noNAN, :);

% Load ADNI2 m12 subinfo
subinfo = csvread([subinfo_dir '/ADNI2_m12_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

% remove ADAS: Recall instruction because CN subjects has 0 std
% remove ANART because it doesn't exist in m6 m12 ...
% combine BNT:without cue and semantic cue
% remove CFT: vegetable because it doesn't exist in ADNI2
% combine MMSE: Language and ADAS: Naming
% combine MMSE: Immediate Recall with MMSE: Delayed Recall
% combine ADAS: Spoken Language, Comprehension
behavior_comb = [behavior(:, 1:4), behavior(:, 5)+2-behavior(:, 19), behavior(:, 6:8), ...
    behavior(:, 10)+behavior(:, 12), behavior(:, [11 13:15]), behavior(:, 16)+behavior(:, 18), ...
    behavior(:, [17 20:21]), sum(behavior(:, 23:24), 2), behavior(:, 27:end)];

% remove behavior with missing scores
ind_noNAN = ~isnan(sum(behavior_comb, 2));
character_rm = character(ind_noNAN, :);
behavior_comb_rm = behavior_comb(ind_noNAN, :);

% remove subjects not from bl
ind_bl_m12 = ismember(character_rm(:, 1), character_rm_rm_ADNI2_bl(:, 1));
character_rm_rm = character_rm(ind_bl_m12, :);
behavior_comb_rm_rm = behavior_comb_rm(ind_bl_m12, :);

% Brain: brain to zscore
ind_brain = CBIG_MMLDA_find_array_in_array(character_rm_rm(:, 1), character(:, 1));
vol4D = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI2_m12/mri/gm_merg_MNI2mm_s10.nii.gz']);
vol4D.vol = vol4D.vol(:, :, :, ind_brain);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz']);
regressor = character(ind_brain, [3 2 5]);
refpara = [out_dir '/ADNI2_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb_rm_rm;
regressor = character_rm_rm(:, [3 2 5]);
refpara = [out_dir '/ADNI2_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

% write zscore to doc
rid = character_rm_rm(:, 1);
ind_out = logical(ones(length(character_rm_rm(:, 4)), 1));
plus_num = 1;
nfold = 1;
out_name = 'ADNI2_m12_ALL_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

%%%
% ADNI1 bl
%%%
subinfo = csvread([subinfo_dir '/ADNI1_bl_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

% remove ADAS: Recall instruction because CN subjects has 0 std
% remove ANART because it doesn't exist in m6 m12 ...
% combine BNT:without cue and semantic cue
% remove CFT: vegetable because it doesn't exist in ADNI2
% combine MMSE: Language and ADAS: Naming
% combine MMSE: Immediate Recall with MMSE: Delayed Recall
% combine ADAS: Spoken Language, Comprehension
behavior_comb = [behavior(:, 1:4), behavior(:, 5)+2-behavior(:, 19), behavior(:, 6:8), ...
    behavior(:, 10)+behavior(:, 12), behavior(:, [11 13:15]), behavior(:, 16)+behavior(:, 18), ...
    behavior(:, [17 20:21]), sum(behavior(:, 23:24), 2), behavior(:, 27:end)];

% remove behavior with missing scores
ind_noNAN = ~isnan(sum(behavior_comb, 2));
character_rm = character(ind_noNAN, :);
behavior_comb_rm = behavior_comb(ind_noNAN, :);

% Brain: create reference paramters 
vol4D_bl = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI1_bl/mri/gm_merg_MNI2mm_s10.nii.gz']);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNI1/ADNI1_bl_gm_mask_MNI2mm.nii.gz']);
ind_CN_bl_brain = (character(:, 4) == 1);
regressor_bl_brain = character(:, [3 2 5]);
zscore_method = 'meanCNstdALL';
output = [out_dir '/ADNI1_bl_brain_refpara.mat'];
CBIG_MMLDA_brain_to_zscore_create_refpara(vol4D_bl, mask, ind_CN_bl_brain, regressor_bl_brain, zscore_method, output);

% Behavior: create reference paramters 
behavior_bl = behavior_comb_rm;
ind_CN_bl_behavior = (character_rm(:, 4) == 1);
regressor_bl_behavior = character_rm(:, [3 2 5]);
zscore_method = 'meanCNstdALL';
output = [out_dir '/ADNI1_bl_behavior_refpara.mat'];
CBIG_MMLDA_behavior_to_zscore_create_refpara(behavior_bl, ind_CN_bl_behavior, ...
    regressor_bl_behavior, zscore_method, output);

% Brain: brain to zscore
ind_brain = CBIG_MMLDA_find_array_in_array(character_rm(:, 1), character(:, 1));
vol4D = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI1_bl/mri/gm_merg_MNI2mm_s10.nii.gz']);
vol4D.vol = vol4D.vol(:, :, :, ind_brain);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNI1/ADNI1_bl_gm_mask_MNI2mm.nii.gz']);
regressor = character(ind_brain, [3 2 5]);
refpara = [out_dir '/ADNI1_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb_rm;
regressor = character_rm(:, [3 2 5]);
refpara = [out_dir '/ADNI1_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

% write zscore to doc
rid = character_rm(:, 1);
ind_out = (character_rm(:, 4) == 3);
plus_num = 1;
nfold = 1;
out_name = 'ADNI1_bl_AD_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

%%%
% ADNI23 PETtau
%%%

% Load ADNI2 bl subinfo to do matrix completion  
subinfo = csvread([subinfo_dir '/ADNI2_bl_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

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
behavior_ADNI2_bl_train = behavior_comb_rm(ind_ALL_noNAN, :);

% Load ADNI23 PETtau subinfo
subinfo = csvread([subinfo_dir '/ADNI23_PETtau_subinfo.csv'], 1, 0);
character = subinfo(:, 1:6);
behavior = subinfo(:, 7:end);

% remove ADAS: Recall instruction because CN subjects has 0 std
% remove ANART because it doesn't exist in m6 m12 ...
% combine BNT:without cue and semantic cue
% remove CFT: vegetable because it doesn't exist in ADNI2
% combine MMSE: Language and ADAS: Naming
% combine MMSE: Immediate Recall with MMSE: Delayed Recall
% combine ADAS: Spoken Language, Comprehension
behavior_comb = [behavior(:, 1:4), behavior(:, 5)+2-behavior(:, 19), behavior(:, 6:8), ...
    behavior(:, 10)+behavior(:, 12), behavior(:, [11 13:15]), behavior(:, 16)+behavior(:, 18), ...
    behavior(:, [17 20:21]), sum(behavior(:, 23:24), 2), behavior(:, 27:end)];

% remove behavior with missing scores > 5
num_of_nan = sum(isnan(behavior_comb), 2);
ind = (num_of_nan <= 5);
character_rm = character(ind, :);
behavior_comb_rm = behavior_comb(ind, :);

% matrix completion for missing scores
behavior_comb_rm_mc = CBIG_MMLDA_matrix_completion_GLM(behavior_comb_rm, behavior_ADNI2_bl_train);

% Brain: brain to zscore
ind_brain = CBIG_MMLDA_find_array_in_array(character_rm(:, 1), character(:, 1));
vol4D = MRIread([ADNI_VBM_path '/preprocessing/output/ADNI23_PETtau/mri/gm_merg_MNI2mm_s10.nii.gz']);
vol4D.vol = vol4D.vol(:, :, :, ind_brain);
mask = MRIread([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz']);
regressor = character(ind_brain, [3 2 5]);
refpara = [out_dir '/ADNI2_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb_rm_mc;
regressor = character_rm(:, [3 2 5]);
refpara = [out_dir '/ADNI2_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

% write zscore to doc
rid = character_rm(:, 1);
ind_out = logical(ones(length(character_rm(:, 4)), 1));
plus_num = 1;
nfold = 1;
out_name = 'ADNI2_PETtau_ALL_meanCNstdALL_plus1';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, out_dir, out_name);

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])