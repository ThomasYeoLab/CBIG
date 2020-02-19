function CBIG_ASDf_FC2doc_estFactors_wrapper(corrMat_ASD_path, corrMat_con_path, sub_info_file, output_dir)
% CBIG_ASDf_FC2doc_estFactors_wrapper(corrMat_ASD_path, corrMat_con_path, sub_info_file, output_dir)
%
% Wrapper function to compute z-normalized, discretized FC and write into
% word documents that will be used by polarLDA model to estimate
% ASD factors. 
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
%     - sub_info_file:
%           Absolute path to the .csv file of subjects' demographic information
%     - output_dir:
%           Absolute path to the output directory where output results will
%           be saved. The outputs are: 
%           1) step1_output_mean_CN.mat: Mean of FC data of control participants
%           2) step1_output_beta_CN.mat: Regression coefficients estimated from control participants
%           3) step1_output_reg_CN_mean_std.mat: Post-regression mean & std of control participants' FC data
%           4) step1_output_zScores.mat: Z-scores & discretized Z-scores of all participants
%           5) step1_output_dx1.dat & step1_output_dx2.dat: "Documents" of
%              ASD and control participants respectively
%
% Example:
%	CBIG_ASDf_FC2doc_estFactors_wrapper('../examples/input/corrMat_ASD_est.mat',
%           '../examples/input/corrMat_Con_est.mat','../examples/input/subInfo_est.csv',
%           '~/Temporary/example_output');
% 
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

%% Concatenate FC matrices together, ASD subjects followed by control subjects
load(corrMat_ASD_path);
load(corrMat_con_path);
corrmat_all = cat(3,corrMat_ASD,corrMat_Con);

%% Get dx_info, cohort_label and regressors
cohort_label = [2; 1]; % 2 represents control, 1 represents ASD

[~, id_dx] = CBIG_ASDf_getSubData(sub_info_file);
dx_info = cell2mat(id_dx(:,2));

regressors = CBIG_ASDf_genRegressors(sub_info_file);

%% Write z-nomralized, discretized FC into word documents
output_name = 'step1_output';
output_dir = [output_dir '/'];
[z, discretized_z] = CBIG_ASDf_FC2doc(corrmat_all, regressors, dx_info, cohort_label, output_dir, output_name);
save([output_dir output_name '_zScores.mat'], 'z', 'discretized_z');

%% Remove path
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
