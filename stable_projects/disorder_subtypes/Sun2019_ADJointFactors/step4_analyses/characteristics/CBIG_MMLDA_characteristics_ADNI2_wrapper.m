function p_all = CBIG_MMLDA_characteristics_wrapper(out_dir)
% p_all = CBIG_MMLDA_characteristics_wrapper(out_dir)
%
% Wrapper function to perform GLM (and logistic regression fro binary variables)
% to compare AD subjects' characteristics (i.e., age, edu, sex, amyloid, apoe2)
% across factors.
%
% Input:
%   - out_dir   : absolute path to the output directory. The output will be a 
%                 csv file like Supplementary Table3 in paper.
%
% Output:
%   - p_all     : all omnibus p values used for FDR correction later
%
% Example:
%   p_all = CBIG_MMLDA_characteristics_wrapper('~/output/analyses/characteristics')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

ADNI_path = [getenv('CBIG_MMLDA_ADNI_DOC_DIR') '/All'];
proj_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];

% get subinfo
subinfo = csvread([proj_dir '/step2_MMLDA/data/ADNI2_bl_subinfo.csv'], 1, 0);

% get factor loadings
rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_AD_meanCNstdALL_plus1_RID.txt'];
inf_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1/' ...
    'k3_inf_ADNI2_bl_AD_meanCNstdALL_plus1-gamma.dat'];
order = [1 3 2];
rid_prob_reorder = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order);

% get age, sex, edu, amyloid, apoe
ind = CBIG_MMLDA_find_array_in_array(rid_prob_reorder(:, 1), subinfo(:, 1));
rid = subinfo(ind, 1);
age = subinfo(ind, 3);
sex = subinfo(ind, 2);
PTDEMOG_file = [ADNI_path '/ADNI_161017/documentation/Subject_Characteristics/PTDEMOG.csv'];
edu = CBIG_MMLDA_get_edu(PTDEMOG_file, repmat({'ADNI2'}, length(rid), 1), ...
    CBIG_MMLDA_matrix2cellstr(rid), repmat({'sc'}, length(rid), 1));
UPENNBIOMK_file = [getenv('CBIG_MMLDA_ADNIMERT_DIR') '/scripts/UPENNBIOMK.csv'];
UCBERKELEYAV45_file = [ADNI_path '/ADNI_180413/documentation/UCBERKELEYAV45_11_14_17.csv'];
[amyloid, ~] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, 'ADNI2', ...
    repmat({'ADNI2'}, length(rid), 1), CBIG_MMLDA_matrix2cellstr(rid), repmat({'bl'}, length(rid), 1));
APOERES_file = [ADNI_path '/ADNI_161017/documentation/Biospecimen/APOERES.csv'];
[apoe2, apoe4] = CBIG_MMLDA_get_apoe(APOERES_file, CBIG_MMLDA_matrix2cellstr(rid));
characteristics = [age edu sex-1 amyloid apoe2 apoe4];

% glm hypothesis test for characteristics
p_all = zeros(size(characteristics, 2), 1);
t = cell(size(characteristics, 2), 5);
for i=1:size(characteristics, 2)
    if i == 3
        [N, wMean, wStdDev, p] = CBIG_MMLDA_characteristics_glm(rid_prob_reorder(:, 2:end), ...
            characteristics(:, i), 'binary');
    else
        [N, wMean, wStdDev, p] = CBIG_MMLDA_characteristics_glm(rid_prob_reorder(:, 2:end), characteristics(:, i));
    end
    p_all(i) = p;
    t{i, 1} = sprintf('N = %d', N);
    t{i, 2} = sprintf('%.3f (%.3f)', wMean(1), wStdDev(1));
    t{i, 3} = sprintf('%.3f (%.3f)', wMean(2), wStdDev(2));
    t{i, 4} = sprintf('%.3f (%.3f)', wMean(3), wStdDev(3));
    t{i, 5} = sprintf('%e', p);
end
disp(t)
cell2csv([out_dir '/characteristics.csv'], t)

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])