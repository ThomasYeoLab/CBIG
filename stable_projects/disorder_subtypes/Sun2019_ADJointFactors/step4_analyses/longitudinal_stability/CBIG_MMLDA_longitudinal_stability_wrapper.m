function CBIG_MMLDA_longitudinal_stability_wrapper(out_dir)
% CBIG_MMLDA_longitudinal_stability_wrapper(out_dir)
%
% Wrapper function to plot longitudinal stability of factors.
%
% Input:
%   - out_dir   : output directory
%
% Example:
%   CBIG_MMLDA_longitudinal_stability_wrapper('~/example')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

%%%
% get factor loadings
%%%
proj_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];
order = [1 3 2];

% AD
rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_AD_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf_ADNI2_bl_AD_meanCNstdALL_plus1-gamma.dat'];
inf1_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf1_ADNI2_bl_AD_meanCNstdALL_plus1-gamma.dat'];
inf2_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf2_ADNI2_bl_AD_meanCNstdALL_plus1-gamma.dat'];
rid_prob_AD_inf = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);
rid_prob_AD_inf1 = CBIG_MMLDA_get_factor_loadings(rid_file, inf1_gamma_file, order);
rid_prob_AD_inf2 = CBIG_MMLDA_get_factor_loadings(rid_file, inf2_gamma_file, order);

% MCI
rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_MCI_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
inf1_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf1_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
inf2_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf2_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
rid_prob_MCI_inf = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);
rid_prob_MCI_inf1 = CBIG_MMLDA_get_factor_loadings(rid_file, inf1_gamma_file, order);
rid_prob_MCI_inf2 = CBIG_MMLDA_get_factor_loadings(rid_file, inf2_gamma_file, order);

% m12 ALL
rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_m12_ALL_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf_ADNI2_m12_ALL_meanCNstdALL_plus1-gamma.dat'];
inf1_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf1_ADNI2_m12_ALL_meanCNstdALL_plus1-gamma.dat'];
inf2_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf2_ADNI2_m12_ALL_meanCNstdALL_plus1-gamma.dat'];
rid_prob_m12_inf = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);
rid_prob_m12_inf1 = CBIG_MMLDA_get_factor_loadings(rid_file, inf1_gamma_file, order);
rid_prob_m12_inf2 = CBIG_MMLDA_get_factor_loadings(rid_file, inf2_gamma_file, order);

% get amyloid
rid_MCI= rid_prob_MCI_inf(:, 1);
ADNI_path = [getenv('CBIG_MMLDA_ANDI_DOC_DIR') '/All'];
UPENNBIOMK_file = [getenv('CBIG_MMLDA_ADNIMERT_DIR') '/scripts/UPENNBIOMK.csv'];
UCBERKELEYAV45_file = [ADNI_path '/ADNI_180413/documentation/UCBERKELEYAV45_11_14_17.csv'];
[~, MCI_state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, 'ADNI2', ...
    repmat({'ADNI2'}, length(rid_MCI), 1), CBIG_MMLDA_matrix2cellstr(rid_MCI), repmat({'bl'}, length(rid_MCI), 1));
rid_Ap_MCI = rid_MCI(MCI_state == 1);
rid_prob_Ap_MCI_inf = rid_prob_MCI_inf(MCI_state==1, :);
rid_prob_Ap_MCI_inf1 = rid_prob_MCI_inf1(MCI_state==1, :);
rid_prob_Ap_MCI_inf2 = rid_prob_MCI_inf2(MCI_state==1, :);

rid_AD = rid_prob_AD_inf(:, 1);
rid_m12 = rid_prob_m12_inf(:, 1);

[rid_AD_inter, ind_bl_AD, ind_m12_AD] = intersect(rid_AD, rid_m12);
[rid_Ap_MCI_inter, ind_bl_Ap_MCI, ind_m12_Ap_MCI] = intersect(rid_Ap_MCI, rid_m12);

% assign colors
colors = [repmat('r', length(rid_AD_inter), 1); repmat('g', length(rid_Ap_MCI_inter), 1)];

% scatter plots
xlabel_name = 'Probability at Baseline';
ylabel_name = 'Probability in Month 12';

prob_x = [rid_prob_AD_inf1(ind_bl_AD, 2:end); rid_prob_Ap_MCI_inf1(ind_bl_Ap_MCI, 2:end)];
prob_y = [rid_prob_m12_inf1(ind_m12_AD, 2:end); rid_prob_m12_inf1(ind_m12_Ap_MCI, 2:end)];
title_name = {'Medial Temporal', 'Lateral Temporal', 'Posterior Cortical'};
out_name = [out_dir '/longitudinal_stability_atrophy'];
CBIG_MMLDA_scatter_plot_k3(prob_x, prob_y, colors, title_name, xlabel_name, ylabel_name, out_name);

prob_x = [rid_prob_AD_inf2(ind_bl_AD, 2:end); rid_prob_Ap_MCI_inf2(ind_bl_Ap_MCI, 2:end)];
prob_y = [rid_prob_m12_inf2(ind_m12_AD, 2:end); rid_prob_m12_inf2(ind_m12_Ap_MCI, 2:end)];
title_name = {'Memory', 'Language', 'Visuospatial Executive Function'};
out_name = [out_dir '/longitudinal_stability_behavior'];
CBIG_MMLDA_scatter_plot_k3(prob_x, prob_y, colors, title_name, xlabel_name, ylabel_name, out_name);

prob_x = [rid_prob_AD_inf(ind_bl_AD, 2:end); rid_prob_Ap_MCI_inf(ind_bl_Ap_MCI, 2:end)];
prob_y = [rid_prob_m12_inf(ind_m12_AD, 2:end); rid_prob_m12_inf(ind_m12_Ap_MCI, 2:end)];
title_name = {'MTL-Memory', 'Lateral Temporal-Language', 'Posterior Cortical Executive'};
out_name = [out_dir '/longitudinal_stability_joint'];
CBIG_MMLDA_scatter_plot_k3(prob_x, prob_y, colors, title_name, xlabel_name, ylabel_name, out_name);

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])