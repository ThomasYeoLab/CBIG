% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% convert vbm and behavior to doc

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

proj_dir = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];
addpath([proj_dir '/step2_MMLDA'])
addpath([proj_dir '/utilities'])

% obtain 28 behavioral description, sign (PLUS:more is worse; MINUS:more is better) and categories in ADNI2
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

% load subinfo
subinfo = csvread(subinfo_path, 1, 0);
character = subinfo(:, 1:3);
behavior = subinfo(:, 4:end);

% load icv
load(gmvol_icv_path)
icv = cell2mat(id_gmVol_icv(:, 3));

% process the behavior
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

% Brain: brain to zscore
vol4D = MRIread(gm_merge_path);
mask = MRIread([proj_dir '/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz']);
regressor = [character(:, [2 3]), icv];
refpara = [proj_dir '/ADNIDataRelease/regression_zscore_paras/ADNI2_bl_brain_refpara.mat'];
Z_brain = CBIG_MMLDA_brain_to_zscore(vol4D,  mask, regressor, refpara);

% Behavior: behavior to zscore
behavior = behavior_comb;
regressor = [character(:, [2 3]), icv];
refpara = [proj_dir '/ADNIDataRelease/regression_zscore_paras/ADNI2_bl_behavior_refpara.mat'];
Z_behavior = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign_comb, refpara);

% write zscore to doc
rid = character(:, 1);
ind_out = logical(ones(length(rid), 1));
plus_num = 1;
nfold = 1;
out_name = 'doc';
CBIG_MMLDA_brain_behavior_zscore_to_doc(Z_brain, Z_behavior, rid, ind_out, plus_num, nfold, doc_out_dir, out_name);

rmpath([proj_dir '/step2_MMLDA'])
rmpath([proj_dir '/utilities'])