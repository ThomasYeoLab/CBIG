% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% IMPORTANT This file is specific to our PNAS paper and our computing cluster. As a result,
% you need to modify the paths at places with comment "please change"

close all; clear; clc;
cd ../lib/

%% Prepare LDA inputs
projDir = '/data/users/xzhang/storage/forPNASRelease/'; % please change
% please change, if you didn't follow this folder structure in VBM
mask = [projDir 'outputs/VBM_bl/concatAndGenMask/GMToNonlinTmp_mod_mean_binThr0.05.nii.gz']; 
inputDir = [getenv('CBIG_CODE_DIR') ...
'/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step2_LDA/replicatePNAS/inputs_brain2doc/'];

% Processing 810 ADNI-1 baseline scans
% please change, if you didn't follow this folder structure in VBM
vol4D = [projDir 'outputs/VBM_bl/smoothing/GMToNonlinTmp_mod_4d_s4.246.nii.gz']; 
% please change, if you didn't follow this folder structure in VBM
concatOrder = [projDir 'outputs/VBM_bl/concatAndGenMask/GMToNonlinTmp_mod_4d_concatOrder.txt']; 
refList = [inputDir 'filenames_228CN.txt'];
nuisanceVars = [inputDir 'age_sex_icv_bl.csv'];
outDir = [projDir 'outputs/LDA_bl/docsForLDA/'];
AD188 = [inputDir 'filenames_188AD.txt'];
MCI394 = [inputDir 'filenames_394MCI.txt'];
CN228 = [inputDir 'filenames_228CN.txt'];
CBIG_brain2doc(vol4D, concatOrder, mask, refList, nuisanceVars, outDir, AD188, MCI394, CN228);

% Processing 560 m24 follow-up scans
% please change, if you didn't follow this folder structure in VBM
vol4D = [projDir 'outputs/VBM_m24/smoothing/GMToNonlinTmp_mod_4d_s4.246.nii.gz']; 
% please change, if you didn't follow this folder structure in VBM
concatOrder = [projDir 'outputs/VBM_m24/concat/GMToNonlinTmp_mod_4d_concatOrder.txt'];
% please change, if you didn't follow this folder structure in VBM
params = [projDir 'outputs/LDA_bl/docsForLDA/refParams.mat']; 
nuisanceVars = [inputDir 'age_sex_icv_m24.csv'];
outDir = [projDir 'outputs/LDA_m24/docsForLDA/'];
m24 = [inputDir 'filenames_560m24.txt']; % since it's the full list, it simply helps reorders the scans
CBIG_brain2doc(vol4D, concatOrder, mask, params, nuisanceVars, outDir, m24);


%% Compile LDA C code
system('rm -rf ./lda-c-dist');
system('rsync -az ${CBIG_CODE_DIR}/external_packages/lda-c-dist ./'); % to keep common space clean
system('cd lda-c-dist; make; cd ..'); % compile the C code


%% Run LDA estimation on 188 AD scans (from 810 ADNI-1 baseline scans)
docs = [projDir 'outputs/LDA_bl/docsForLDA/docs_AD188.dat'];
R = 20; % number of random initializations
outDir = [projDir 'outputs/LDA_bl/est/'];
% (optional) job queue; If you have a cluster, use it to specify a job queue, initializations will run in parallel
queue = 'circ-spool'; 
% Try K = 2, 3 and 4
for K = 2:4
    if isempty(queue)
        cmd = ['./CBIG_LDA_est.sh -d ' docs ' -k ' num2str(K) ' -r ' num2str(R) ' -o ' outDir];
    else
        cmd = ['./CBIG_LDA_est.sh -d ' docs ' -k ' num2str(K) ' -r ' num2str(R) ' -o ' outDir ' -q ' queue];
    end
    system(cmd);
end
% Hold off the next section untill all jobs are finished
CBIG_waitUntilFinished(progressFile, noJobs);


%% Visualize the estimated topics
GMMask = [projDir 'outputs/VBM_bl/concatAndGenMask/GMToNonlinTmp_mod_mean_binThr0.05.nii.gz'];
for K = 2:4
    initFolders = [projDir 'outputs/LDA_bl/est/K' num2str(K) '/R*'];
    outDir = [projDir 'outputs/LDA_bl/results/K' num2str(K) '/'];
    if K == 3 % save model name for future inference
        model = CBIG_visualizeFactors(initFolders, GMMask, outDir);
    else % discard model names
        CBIG_visualizeFactors(initFolders, GMMask, outDir);
    end
end


%% Run LDA inference on 394 MCI, 228 CN and 560 m24 follow-up scans
% 394 MCI
docs = [projDir 'outputs/LDA_bl/docsForLDA/docs_MCI394.dat'];
outName = [projDir 'outputs/LDA_bl/inf/K3_394MCI'];
cmd = ['./CBIG_LDA_inf.sh -d ' docs ' -m ' model ' -o ' outName];
system(cmd);
% 228 CN
docs = [projDir 'outputs/LDA_bl/docsForLDA/docs_CN228.dat'];
outName = [projDir 'outputs/LDA_bl/inf/K3_228CN'];
cmd = ['./CBIG_LDA_inf.sh -d ' docs ' -m ' model ' -o ' outName];
system(cmd);
% 560 m24
docs = [projDir 'outputs/LDA_m24/docsForLDA/docs_m24.dat'];
outName = [projDir 'outputs/LDA_m24/inf/K3_m24'];
cmd = ['./CBIG_LDA_inf.sh -d ' docs ' -m ' model ' -o ' outName];
system(cmd);
