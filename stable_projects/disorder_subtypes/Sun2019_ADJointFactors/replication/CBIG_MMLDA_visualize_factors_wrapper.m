function CBIG_MMLDA_visualize_factors_wrapper(mmlda_dir, visualize_dir)
% This wrapper visualize factors after MMLDA estimation
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
unit_test_path = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];

if nargin < 1
    mmlda_dir = [unit_test_path '/step2_MMLDA/results/estimation'];
end
if nargin < 2
    visualize_dir = [unit_test_path '/step2_MMLDA/results/visualizeFactors'];
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])

t = readtable([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors' ...
    '/utilities/Behavior_Abbreviation.csv']);
behavior_name = t.Behavioral_Name_Short;
behavior_domain = t.Behavioral_Cate;

% ADNI2 bl K=2,3,4
in_dir = [mmlda_dir '/ADNI2_bl_AD_meanCNstdALL_plus1'];
out_dir = [visualize_dir '/ADNI2_bl_AD_meanCNstdALL_plus1'];
mask = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz'];
min_thresh = 7.5e-6;
max_thresh = 1.5e-5;
for k = 2:4
    CBIG_MMLDA_visualize_factors(in_dir, out_dir, k, min_thresh, max_thresh, mask, behavior_name, behavior_domain)
end

% ADNI2 bl K=3 10folds
for i = 1:10
    in_dir = [mmlda_dir '/ADNI2_bl_AD_meanCNstdALL_plus1_train' num2str(i)];
    out_dir = [visualize_dir '/ADNI2_bl_AD_meanCNstdALL_plus1_train' num2str(i)];
    mask = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNIGO2/ADNI2_bl_gm_mask_MNI2mm.nii.gz'];
    min_thresh = 7.5e-6;
    max_thresh = 1.5e-5;
    k = 3;
    CBIG_MMLDA_visualize_factors(in_dir, out_dir, k, min_thresh, max_thresh, mask, behavior_name, behavior_domain)
end

% ADNI1 bl K=3
in_dir = [mmlda_dir '/ADNI1_bl_AD_meanCNstdALL_plus1'];
out_dir = [visualize_dir '/ADNI1_bl_AD_meanCNstdALL_plus1'];
mask = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
    'Sun2019_ADJointFactors/ADNIDataRelease/SPM_VBM_files/ADNI1/ADNI1_bl_gm_mask_MNI2mm.nii.gz'];
min_thresh = 7.5e-6;
max_thresh = 1.5e-5;
k = 3;
CBIG_MMLDA_visualize_factors(in_dir, out_dir, k, min_thresh, max_thresh, mask, behavior_name, behavior_domain)

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/step2_MMLDA'])

