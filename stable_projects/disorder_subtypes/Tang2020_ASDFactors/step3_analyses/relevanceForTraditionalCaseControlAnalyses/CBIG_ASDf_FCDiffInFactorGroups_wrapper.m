function p_NBS = CBIG_ASDf_FCDiffInFactorGroups_wrapper(outputDir, cluster)
% p_NBS = CBIG_ASDf_FCDiffInFactorGroups_wrapper(outputDir, cluster)
% 
% Wrapper function to compare FC difference between all ASD and control subjects,
% as well as between ASD and control subjects within each factor group on
% ABIDE-I sample. Statistical significance is tested using network-based
% statistic (NBS; Zalesky et al., 2010). FC differences will also be
% plotted as 419x419 matrices (NBS thresholded).
% In our paper, we experimented with two different  subgrouping criteria:
% (A) If an ASD subject's Pr(Factor|Participant) > 0.5, he/she was assigned to
% the corresponding factor group;
% (B) The ASD subjects were assigned to the corresponding factor with highest
% Pr(Factor|Participant).
% 
% NOTE: NBS takes some time to run. Alternatively, if you have access to
% our circ-spool cluster, you can submit jobs to run in parallel.
%
% Input:
%     - outputDir:
%           Absolute path to the directory where output will be saved
%     - cluster:
%           String. Cluster name. If no cluster, enter [] or ''.
% Output:
%      - p_NBS:
%           All p-values from network-based statistic analysis. These
%           p-values can be used later for FDR multiple comparisons
%           correction.
% 
% Example:
%       p_NBS = CBIG_ASDf_FCDiffInFactorGroups_wrapper('~/example_output/step3_analyses', 
%       'circ-spool')
%       Will run the NBS by submitting jobs to 'circ-spool'
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

Nperm = 10000;
tThresh = 2;  % t-statistic threshold used for plots in our paper
p_NBS = [];

%% Define and add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
addpath(fullfile(CODE_DIR,'step2_polarLDA'));

corrMat_dir = fulfile('mnt','eql','yeo11','data','ABIDE1_preprocess','fc_correlation','output_316');
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
inputDir = fullfile(UNIT_TEST_DIR,'data');

sub_info_file = fullfile(inputDir,'subInfo_316.csv');
subgrpDir = fullfile(inputDir,'ABIDE1_subgroups');

%% Load ABIDE-I ASD participants' IDs in each sub-group
id_asd_file = fullfile(subgrpDir,'id_asd.mat');
id_con_file = fullfile(subgrpDir,'id_con.mat');

id_K3F1_file = fullfile(subgrpDir,'id_K3F1.mat');
id_K3F2_file = fullfile(subgrpDir,'id_K3F2.mat');
id_K3F3_file = fullfile(subgrpDir,'id_K3F3.mat');
id_K3Con1_file = fullfile(subgrpDir,'id_K3Con1.mat');
id_K3Con2_file = fullfile(subgrpDir,'id_K3Con2.mat');
id_K3Con3_file = fullfile(subgrpDir,'id_K3Con3.mat');

%% Set random number seed
rng('default');

%% FC difference between all ASD and all control subjects
disp('----FC difference between all ASD and all control subjects:');
output_name = ['nbs_allSub_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];

%%% Run NBS, either submit job to CBIG cluster or run directly on matlab
if ~isequal(cluster,'circ-spool')
    curr_output_name = fullfile(outputDir, output_name);
    p_all = CBIG_ASDf_FCDiffAllSub_NBS(tThresh, id_asd_file, id_con_file, sub_info_file, ...
Nperm, curr_output_name);
    p_NBS = [p_NBS; p_all];
else
    cmd = ['./CBIG_ASDf_NBSall_job.sh ' num2str(tThresh) ' ' id_asd_file ' ' ...
id_con_file ' ' sub_info_file ' ' num2str(Nperm) ' ' output_name ' ' outputDir];
    system(cmd);
end

%% FC difference in factor groups, K = 3
disp('----FC difference between ASD and control subjects within factor group (K = 3):');

%%% K = 3, factor 1 subgroup
output_name = ['nbs_K3F1_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
if ~isequal(cluster,'circ-spool')
    curr_output_name = fullfile(outputDir, output_name);
    p_F1 = CBIG_ASDf_FCDiffInSubgrp_NBS(tThresh, id_K3F1_file, id_K3Con1_file, id_asd_file, ...
id_con_file, sub_info_file, Nperm, curr_output_name);
    p_NBS = [p_NBS; p_F1];
else
    cmd = ['./CBIG_ASDf_NBSsubgrp_job.sh ' num2str(tThresh) ' ' id_K3F1_file ' ' ...
id_K3Con1_file ' ' id_asd_file ' ' id_con_file ' ' sub_info_file ' ' num2str(Nperm) ' ' ...
output_name ' ' outputDir];
    system(cmd);
end

%%% K = 3, factor 2 subgroup
output_name = ['nbs_K3F2_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
if ~isequal(cluster,'circ-spool')
    curr_output_name = fullfile(outputDir, output_name);
    p_F2 = CBIG_ASDf_FCDiffInSubgrp_NBS(tThresh, id_K3F2_file, id_K3Con2_file, id_asd_file, ...
id_con_file, sub_info_file, Nperm, curr_output_name);
    p_NBS = [p_NBS; p_F2];
else
    cmd = ['./CBIG_ASDf_NBSsubgrp_job.sh ' num2str(tThresh) ' ' id_K3F2_file ' ' ...
id_K3Con2_file ' ' id_asd_file ' ' id_con_file ' ' sub_info_file ' ' num2str(Nperm) ' ' ...
output_name ' ' outputDir];
    system(cmd);
end

%%% K = 3, factor 3 subgroup
output_name = ['nbs_K3F3_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
if ~isequal(cluster,'circ-spool')
    curr_output_name = fullfile(outputDir, output_name);
    p_F3 = CBIG_ASDf_FCDiffInSubgrp_NBS(tThresh, id_K3F3_file, id_K3Con3_file, id_asd_file, ...
id_con_file, sub_info_file, Nperm, curr_output_name);
    p_NBS = [p_NBS; p_F3];
else
    cmd = ['./CBIG_ASDf_NBSsubgrp_job.sh ' num2str(tThresh) ' ' id_K3F3_file ' ' ...
id_K3Con3_file ' ' id_asd_file ' ' id_con_file ' ' sub_info_file ' ' num2str(Nperm) ' ' ...
output_name ' ' outputDir];
    system(cmd);
end

%% Plot FC differences, thresholded by NBS
%%% All ASD vs controls
if ~isempty(cluster) % If submitted jobs, wait until jobs finished
    CBIG_ASDf_checkJobStatus(fullfile(outputDir,'progressFile.txt'), 4, 600);
end
output_name = ['nbs_allSub_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
load(fullfile(outputDir, [output_name '.mat'])) % load NBS result
CBIG_ASDf_plotFCDiff_NBSThresholded(sub_info_file, id_asd_file, id_con_file, ADJ, 1, ...
corrMat_dir, [-0.08 0.08], fullfile(outputDir, output_name));

%%% K = 3, factor 1 subgroup
output_name = ['nbs_K3F1_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
load(fullfile(outputDir, [output_name '.mat'])) % load NBS result
CBIG_ASDf_plotFCDiff_NBSThresholded(sub_info_file, id_K3F1_file, id_K3Con1_file, ADJ, 1, ...
corrMat_dir, [-0.08 0.08], fullfile(outputDir, output_name));

%%% K = 3, factor 2 subgroup
output_name = ['nbs_K3F2_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
load(fullfile(outputDir, [output_name '.mat'])) % load NBS result
CBIG_ASDf_plotFCDiff_NBSThresholded(sub_info_file, id_K3F2_file, id_K3Con2_file, ADJ, 1, ...
corrMat_dir, [-0.08 0.08], fullfile(outputDir, output_name));


%%% K = 3, factor 3 subgroup
output_name = ['nbs_K3F3_Nperm' num2str(Nperm) '_tThresh' num2str(tThresh)];
load(fullfile(outputDir, [output_name '.mat'])) % load NBS result
CBIG_ASDf_plotFCDiff_NBSThresholded(sub_info_file, id_K3F3_file, id_K3Con3_file, ADJ, 1, ...
corrMat_dir, [-0.08 0.08], fullfile(outputDir, output_name));

%% Remove paths
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
rmpath(fullfile(CODE_DIR,'step2_polarLDA'));
