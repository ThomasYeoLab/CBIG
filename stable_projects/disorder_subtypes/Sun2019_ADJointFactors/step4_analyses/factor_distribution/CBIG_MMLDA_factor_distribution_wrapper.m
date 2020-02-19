function CBIG_MMLDA_factor_distribution_wrapper(out_dir)
% CBIG_MMLDA_factor_distribution_wrapper(out_dir)
%
% Wrapper function to plot factor distribution of ADNIGO/2 AD and A+MCI subjects.
%
% Input:
%   - out_dir   : output directory
%
% Example:
%   CBIG_MMLDA_factor_distribution_wrapper('~/example')
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
rid_prob_AD_inf = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_AD_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf1_ADNI2_bl_AD_meanCNstdALL_plus1-gamma.dat'];
rid_prob_AD_inf1 = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_AD_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf2_ADNI2_bl_AD_meanCNstdALL_plus1-gamma.dat'];
rid_prob_AD_inf2 = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

% MCI
rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_MCI_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
rid_prob_MCI_inf = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_MCI_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf1_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
rid_prob_MCI_inf1 = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_MCI_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf2_ADNI2_bl_MCI_meanCNstdALL_plus1-gamma.dat'];
rid_prob_MCI_inf2 = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

rid_AD = rid_prob_AD_inf(:, 1);
rid_MCI= rid_prob_MCI_inf(:, 1);
ADNI_path = [getenv('CBIG_MMLDA_ADNI_DIR') '/All'];
UPENNBIOMK_file = [getenv('CBIG_MMLDA_ADNIMERT_DIR') '/scripts/UPENNBIOMK.csv'];
UCBERKELEYAV45_file = [ADNI_path '/ADNI_180413/documentation/UCBERKELEYAV45_11_14_17.csv'];
[~, AD_state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, 'ADNI2', ...
    repmat({'ADNI2'}, length(rid_AD), 1), CBIG_MMLDA_matrix2cellstr(rid_AD), repmat({'bl'}, length(rid_AD), 1));
[~, MCI_state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, 'ADNI2', ...
    repmat({'ADNI2'}, length(rid_MCI), 1), CBIG_MMLDA_matrix2cellstr(rid_MCI), repmat({'bl'}, length(rid_MCI), 1));
rid_Ap_MCI = rid_MCI(MCI_state == 1);
rid_prob_Ap_MCI_inf = rid_prob_MCI_inf(MCI_state==1, :);
rid_prob_Ap_MCI_inf1 = rid_prob_MCI_inf1(MCI_state==1, :);
rid_prob_Ap_MCI_inf2 = rid_prob_MCI_inf2(MCI_state==1, :);

% AD triangle plot
% joint
colors = [];
for i = 1:size(rid_prob_AD_inf, 1)
    if AD_state(i) == 1
        colors = [colors; 'r'];
    elseif AD_state(i) == 0
        colors = [colors; 'g'];
    elseif isnan(AD_state(i))
        colors = [colors; 'b'];
    end
end
vertex_label = {'MT-M', 'LT-L', 'PC-E'};
bary_coor = rid_prob_AD_inf(:, 2:end);
out_name = [out_dir '/AD_factor_distribution_joint'];
CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name);

% atrophy
bary_coor = rid_prob_AD_inf1(:, 2:end);
out_name = [out_dir '/AD_factor_distribution_atrophy'];
CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name);

% behavior
bary_coor = rid_prob_AD_inf2(:, 2:end);
out_name = [out_dir '/AD_factor_distribution_behavior'];
CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name);

% A+MCI triangle plot
% joint
bary_coor = rid_prob_Ap_MCI_inf(:, 2:end);
colors = [repmat('r', size(bary_coor, 1), 1)];
vertex_label = {'MT-M', 'LT-L', 'PC-E'};
out_name = [out_dir '/A+MCI_factor_distribution_joint'];
CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name);

% atrophy
bary_coor = rid_prob_Ap_MCI_inf1(:, 2:end);
colors = [repmat('r', size(bary_coor, 1), 1)];
vertex_label = {'MT-M', 'LT-L', 'PC-E'};
out_name = [out_dir '/A+MCI_factor_distribution_atrophy'];
CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name);

% behavior
bary_coor = rid_prob_Ap_MCI_inf2(:, 2:end);
colors = [repmat('r', size(bary_coor, 1), 1)];
vertex_label = {'MT-M', 'LT-L', 'PC-E'};
out_name = [out_dir '/A+MCI_factor_distribution_behavior'];
CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name);

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])