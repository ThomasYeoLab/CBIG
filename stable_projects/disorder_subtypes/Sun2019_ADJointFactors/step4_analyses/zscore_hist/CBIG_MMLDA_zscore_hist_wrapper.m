function CBIG_MMLDA_zscore_hist_wrapper(out_dir)
% CBIG_MMLDA_zscore_hist_wrapper(out_dir)
%
% Wrapper function to plot the histogram of zscores.
%
% Input:
%   - out_dir   : output directory
%
% Example:
%   CBIG_MMLDA_zscore_hist_wrapper('~/example')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])

subinfo_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/step2_MMLDA/data'];
ADNI_VBM_path = [getenv('CBIG_MMLDA_ANDI_DIR') '/Sun2019_SPMVBM'];

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

%%%
% hist of zscores
%%%
bin_vector = [-5.5:1:7.5];

figure
% figure size
set(gcf,'Position',[100, 600, 600, 450]);

% bar face and edge color
[N, X] = hist(Z_behavior(:), bin_vector);
[~, ind] = max(N);
h = bar(X, N, 'hist');
set(h, 'facecolor', 'r', 'edgecolor', 'w')

% figure axis
xlim([-5 7])
set(gca, 'xtick', [-5:1:7], 'fontsize', 20)
set(gca, 'linewidth', 1)
axis square;
box off;

% figure text
xlabel('Z score')
ylabel('Count')
title('ADNI2 bl behavior')

hgexport(gcf, [out_dir '/ADNIGO2_bl_behavior'])
eps2xxx([out_dir '/ADNIGO2_bl_behavior.eps'], {'png'})


figure
% figure size
set(gcf,'Position',[100, 600, 600, 450]);

% bar face and edge color
[N, X] = hist(Z_brain(:), bin_vector);
[~, ind] = max(N);
h = bar(X, N, 'hist');
set(h, 'facecolor', 'r', 'edgecolor', 'w')

% figure axis
xlim([-5 7])
set(gca, 'xtick', [-5:1:7], 'fontsize', 20)
set(gca, 'linewidth', 1)
axis square;
box off;

% figure text
xlabel('Z score')
ylabel('Count')
title('ADNI2 bl atrophy')

hgexport(gcf, [out_dir '/ADNIGO2_bl_atrophy'])
eps2xxx([out_dir '/ADNIGO2_bl_atrophy.eps'], {'png'})

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])