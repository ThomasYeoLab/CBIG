function CBIG_ASDf_FC2doc_infFactorComp_wrapper(corrMat_ASD_path, corrMat_con_path, ref_path, sub_info_file, output_dir)
% CBIG_ASDf_FC2doc_infFactorComp_wrapper(corrMat_ASD_path, corrMat_con_path, ref_path, sub_info_file, output_dir)
%
% Wrapper function to compute Z-normalized, discretized FC (Z-normalization
% w.r.t controls in ABIDE-II+GENDAAR) and write into word documents 
% that will be used by polarLDA model to infer factor compositions 
% of ASD subjects in ABIDE-I.
%
% Input:
%     - corrMat_ASD_path:
%           Absolute path to the functional connectivity matrices (.mat file) of ASD
%           subjects. Assuming .mat file is named as "corrMat_ASD", and it is
%           of size MxMxN1, where M is the number of ROIs, N1 is the number
%           of ASD participants.
%     - corrMat_con_path:
%           Absolute path to the functional connectivity matrices (.mat file) of 
%           control subjects. Assuming .mat file is named as "corrMat_Con", and it is
%           of size MxMxN2, where M is the number of ROIs, N2 is the number
%           of control participants.
%     - ref_path:
%           Absolute path to the mean and std of controls' FC data in
%           ABIDE-II+GENDAAR. i.e., the absolute path to the output file
%           output_name_reg_CN_mean_std.mat from function
%           CBIG_ASDf_FC2doc.m
%     - sub_info_file:
%           Absolute path to the .csv file of subjects' demographic information
%     - output_dir:
%           Absolute path to the output directory where output results will
%           be saved. The outputs are: 
%           1) step1_output_inf_mean_CN.mat: Mean of FC data of control participants
%           2) step1_output_inf_beta_CN.mat: Regression coefficients estimated from control participants
%           4) step1_output_inf_zScores.mat: Z-scores & discretized Z-scores of all participants
%           5) step1_output_inf_dx1.dat & step1_output_inf_dx2.dat: "Documents" 
%              of ASD and control participants respectively
%
% Example:
%	CBIG_ASDf_FC2doc_infFactorComp_wrapper('../examples/input/corrMat_ASD_inf.mat',
%           '../examples/input/corrMat_Con_inf.mat',
%           '~/Temporary/example_output/step1_output_reg_CN_mean_std.mat',
%           '../examples/input/subInfo_inf.csv','~/Temporary/example_output')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

%% Load FC matrices and concatenate together, ASD followed by control subjects
load(corrMat_ASD_path);
load(corrMat_con_path);
corrMat_all = cat(3,corrMat_ASD,corrMat_Con);

%% Get dx_info, cohort_label and regressors
cohort_label = [2; 1]; % 2 represents control, 1 represents ASD

[~, id_dx] = CBIG_ASDf_getSubData(sub_info_file);
dx_info = cell2mat(id_dx(:,2));

regressors = CBIG_ASDf_genRegressors(sub_info_file);

%% Load mean and std of post-regression FC values of control subjects in ABIDE-II+GENDAAR
load(ref_path);

%% Write correlation into word documents
output_name = 'step1_output_inf';
output_dir = [output_dir '/'];
[z, discretized_z] = CBIG_ASDf_FC2doc_forInference(corrMat_all, ...
reg_CN_mean, reg_CN_std, regressors, dx_info, cohort_label, output_dir, output_name);
save(fullfile(output_dir, [output_name '_zScores.mat']), 'z', 'discretized_z');

%% Remove path
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
