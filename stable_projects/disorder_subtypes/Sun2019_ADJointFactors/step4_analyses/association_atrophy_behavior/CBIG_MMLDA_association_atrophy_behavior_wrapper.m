function p_all = CBIG_MMLDA_association_atrophy_behavior_wrapper(out_dir)
% CBIG_MMLDA_association_atrophy_behavior_wrapper(out_dir)
%
% Wrapper function to use GLM to investigate association between atrophy and behavior
%
% Input:
%   - out_dir   : output directory
%
% Output:
%   - p_all     : all pairwise p values
%
% Example:
%   CBIG_MMLDA_association_atrophy_behavior_wrapper('~/example')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(genpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities']))

%%%
% get factor loadings
%%%
proj_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];
order = [1 3 2];

rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_MCI_meanCNstdALL_plus1_RID.txt'];
inf1_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf1_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
inf2_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf2_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
rid_prob_MCI_inf1 = CBIG_MMLDA_get_factor_loadings(rid_file, inf1_gamma_file, order);
rid_prob_MCI_inf2 = CBIG_MMLDA_get_factor_loadings(rid_file, inf2_gamma_file, order);

rid_MCI= rid_prob_MCI_inf1(:, 1);
ADNI_path = [getenv('CBIG_MMLDA_ANDI_DOC_DIR') '/All'];
UPENNBIOMK_file = [getenv('CBIG_MMLDA_ADNIMERT_DIR') '/scripts/UPENNBIOMK.csv'];
UCBERKELEYAV45_file = [ADNI_path '/ADNI_180413/documentation/UCBERKELEYAV45_11_14_17.csv'];
[~, MCI_state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, 'ADNI2', ...
    repmat({'ADNI2'}, length(rid_MCI), 1), CBIG_MMLDA_matrix2cellstr(rid_MCI), repmat({'bl'}, length(rid_MCI), 1));
rid_Ap_MCI = rid_MCI(MCI_state == 1);
rid_prob_Ap_MCI_inf1 = rid_prob_MCI_inf1(MCI_state==1, :);
rid_prob_Ap_MCI_inf2 = rid_prob_MCI_inf2(MCI_state==1, :);

%%%
% get age sex
%%%
subinfo = csvread([proj_dir '/step2_MMLDA/data/ADNI2_bl_subinfo.csv'], 1, 0);
ind = CBIG_MMLDA_find_array_in_array(rid_Ap_MCI, subinfo(:, 1));
rid = subinfo(ind, 1);
age = subinfo(ind, 3);
sex = subinfo(ind, 2);

%%%
% GLM fit, hypothese test and forest plot 
%%%
k = 3;
X = [rid_prob_Ap_MCI_inf1(:, 3:end), age, sex];
dependency = 'dependent';
compare_factor_names = {'L-M', 'C-M', 'C-L'};
CI_scale = 1;
p_thresh = 0.0167;

y = rid_prob_Ap_MCI_inf2(:, 2);
out_name = 'association_atrophy_memory_A+MCI';
p_memory = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, ...
    out_name, compare_factor_names, CI_scale, p_thresh);

y = rid_prob_Ap_MCI_inf2(:, 3);
out_name = 'association_atrophy_language_A+MCI';
p_language = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, ...
    out_name, compare_factor_names, CI_scale, p_thresh);

y = rid_prob_Ap_MCI_inf2(:, 4);
out_name = 'association_atrophy_ef_A+MCI';
p_ef = CBIG_MMLDA_glmfit_hypotest_forestplot(k, X, y, dependency, out_dir, ...
    out_name, compare_factor_names, CI_scale, p_thresh);

p_all = [p_memory; p_language; p_ef];

rmpath(genpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities']))