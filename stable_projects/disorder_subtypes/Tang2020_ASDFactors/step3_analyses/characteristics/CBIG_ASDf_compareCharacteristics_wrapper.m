function p_characteristics = CBIG_ASDf_compareCharacteristics_wrapper(factorLoading_dir, outputDir)
% p_characteristics = CBIG_ASDf_compareCharacteristics_wrapper(factorLoading_dir, outputDir)
% 
% Wrapper function to perform GLM (and logistic regression for binary
% variables) to compare ASD subjects' characteristics (i.e., age,
% sex, IQ and head motion) across latent factors.
% 
% Input:
%     - factorLoading_dir:
%           Absolute path to "factorComp.txt" file, where each row
%           corresponds to the factor loadings of one ASD participant.
%     - outputDir:
%           Absolute path to the output directory
% Output:
%     - p_characteristics:
%           All pairwise p-values obtained from all GLM and logistic
%           regressions. These p-values can be used for FDR multiple
%           comparisons correction later.
%
% Example:
%       p = CBIG_ASDf_compareCharacteristics_wrapper('~/example_output/step2_polarLDA/factorComp.txt',
%       '~/example_output/analyses/characteristics')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

p_characteristics = [];

%% Add path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','characteristics'));
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
inputDir = fullfile(UNIT_TEST_DIR,'data');
factorOrder = [3 2 1];
k = 3;

sub_info_file = fullfile(inputDir,'subInfo_654.csv');
load(fullfile(inputDir,'id_sexBalanced.mat')); % Load subject IDs for sex difference analysis (114 subjects)

%% Get diagnosis, age, motion, FIQ scores
[~, id_dx, id_age, ~, id_motion, id_fiq, ~, ~] = CBIG_ASDf_getSubData(sub_info_file);

dx_info = cell2mat(id_dx(:,2));
age = cell2mat(id_age(:,2));
motion = cell2mat(id_motion(:,2));
FIQ_score = cell2mat(id_fiq(:,2));

%% Get regressors (i.e., age, sex, motion and sites)
regressors = CBIG_ASDf_genRegressors(sub_info_file);
idx_asd = (dx_info == 1);
reg = regressors(idx_asd,4:end); % Sites as regressors

%% Construct GLMs 
factorComp = dlmread(factorLoading_dir);
factorComp = factorComp(:,factorOrder);
if size(factorComp,2) ~= 3
    error('NOT CONFIGURED! Provide factor loading file for K=3.');
end

disp('----------------K = 3:');
%%% Compare age across factors
disp('----AGE:');
X = [factorComp(:,1:(k-1)) reg];
y = age(idx_asd,1);
numReg = size(reg,2);
output_name = fullfile(outputDir, ['k' num2str(k) '_AGE']);
curr_p = CBIG_ASDf_fitGLM_hypoTest(k, X, y, numReg, output_name); 
p_characteristics = [p_characteristics; curr_p(2:end,1)];

%%% Compare motion across factors
disp('----MOTION:');
X = [factorComp(:,1:(k-1)) reg];
y = motion(idx_asd,1);
numReg = size(reg,2);
output_name = fullfile(outputDir, ['k' num2str(k) '_MOTION']);
curr_p = CBIG_ASDf_fitGLM_hypoTest(k, X, y, numReg, output_name);
p_characteristics = [p_characteristics; curr_p(2:end,1)];

%%% Compare FIQ across factors
disp('----FIQ:');
X = [factorComp(:,1:(k-1)) reg];
y = FIQ_score(idx_asd,1);
numReg = size(reg,2);
output_name = fullfile(outputDir, ['k' num2str(k) '_FIQ']);
curr_p = CBIG_ASDf_fitGLM_hypoTest(k, X, y, numReg, output_name);
p_characteristics = [p_characteristics; curr_p(2:end,1)];

%%% Compare sex across factors
disp('----SEX:');
output_name = fullfile(outputDir, 'k' num2str(k) '_SEX']);
[~, curr_out] = CBIG_ASDf_logReg_compareSex(k, id_sexBalanced, ...
sub_info_file, factorLoading_dir, factorOrder, output_name);
p_characteristics = [p_characteristics; curr_out(2:end,1)];

%% Remove path
rmpath(fullfile(CODE_DIR,'step3_analyses','characteristics'));
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

