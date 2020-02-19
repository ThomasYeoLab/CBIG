function CBIG_VK2019_replication_wrapper
%
% This wrapper runs the main analyses performed in our Biological Psychiatry
% paper. The steps are:
% 
% 1. Loading RSFC & behavior data
% 2. Regressing out confounds from RSFC & behavior data
% 3. PLS analysis
% 4. Permutation testing over PLS singular values
% 5. Plotting RSFC & behavior composite scores & loadings 
% 6. Bootstrap resampling to test significance of RSFC & behavior loadings
% 7. PLS in cross-validation
% 8. Validation of latent components using task FC data
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

current_dir = fileparts(mfilename('fullpath'));
pos_v = strfind(current_dir,filesep);
root_dir = fullfile(current_dir(1:pos_v(length(pos_v)) - 1));
scripts_dir = fullfile(root_dir,'replication','code');
out_dir = fullfile(root_dir,'replication','output');
if ~exist(out_dir), mkdir(out_dir); end

addpath(scripts_dir);
addpath(fullfile(getenv('CBIG_CODE_DIR'),'external_packages','matlab','non_default_packages','PLS_MIPlab'));

data_root_dir = getenv('CBIG_VK2019_UCLA_CNP_DIR');
behav_dir = fullfile(data_root_dir,'behavData','phenotype','mat');
motion_dir = fullfile(data_root_dir,'preprocessedData','rsfMRI','GSR');
FC_dir = fullfile(data_root_dir,'FC','rest','fsavg6_censor_GSR');
taskFC_dir = fullfile(data_root_dir,'FC','averageAcrossTasks','fsavg6_censor_GSR');

subjListFile = fullfile(root_dir,'data_release','PLS_loadings','LC_composite_scores.csv');
subjListFile_all = fullfile(data_root_dir,'behavData','diagnosis.csv');


% Options
nPerms_rest = 1000; % permutations in main (RSFC) PLS analysis
nPerms_cv = 1000; % permutations in cross-validation
nPerms_task = 1000; % permutations in task FC validation
nBootstraps = 500; % bootstraps for RSFC & behavior loadings
normalization_img = 1; % normalization options for RSFC data
normalization_behav = 1; % normalization options for behavior adta
% 0 = no normalization
% 1 = zscore across all subjects
% 2 = zscore within groups (default)
% 3 = std normalization across subjects (no centering)
% 4 = std normalization within groups (no centering)

confounds= {'age','gender','educ1','scanner','motion1'}; % Regressors (motion1=FDRMS; motion2=FDRMS & DVARS)
RSfile_stem = 'bld001_rest_skip4_stc';

%% 1. Load data

disp('(1) Loading data');

% Load list of subjects & diagnoses
fid = fopen(subjListFile_all);
file = textscan(fid,'%s%s','delimiter',',');
fclose(fid);
subj_names = file{1};
clear file fid

% Load list of included subjects
fid = fopen(subjListFile);
file = textscan(fid,'%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
namesInclSubj = file{1}(2:end,:);
nSubj = numel(namesInclSubj);
diagInclSubj = file{2}(2:end);
CONST_DIAGNOSIS = unique(diagInclSubj);
clear file fid

% Get index of included subjects & code diagnosis as numerical value
commonInclSubj = nan(nSubj,1);
diagnosis_grouping = nan(nSubj,1);

for iter_subj = 1:nSubj
    a = strfind(subj_names,namesInclSubj{iter_subj});
    commonInclSubj(iter_subj,1) = find(not(cellfun('isempty',a)));     
    b = strcmp(CONST_DIAGNOSIS,diagInclSubj{iter_subj});
    diagnosis_grouping(iter_subj,1) = find(b==1);     
    clear a b
end

% Behavior variables
behavTests = {'asrs','hopkins','barratt','dickman','mpq','eysenck','bipolarII_scales','golden',...
    'hypomanic','chapman','tci','cvlt','wms','wais','color_trails','dkefs_english','task_switch','ant','cpt'};

% Load behavior scores
this_num = 1;
for iter_behav = 1:numel(behavTests)
    foo = load(fullfile(behav_dir,behavTests{iter_behav})); 
    
    for iter_score = 1:size(foo.scoreNames,2)      
        behavData(:,this_num) = foo.score(:,iter_score);
        behavNames{this_num,1} = foo.scoreNames{1,iter_score};
        this_num = this_num + 1;
    end
    clear foo
end

% Store behavior data of all included subjects in Y
for iter_behav = 1:numel(behavNames)
    for iter_subj = 1:nSubj
        Y_nonreg(iter_subj,iter_behav) = behavData(commonInclSubj(iter_subj),iter_behav);
    end
end

% Load FC matrices, vectorize them and put them in X
for iter_subj = 1:nSubj
    clear RSfile foo
    RSfile = dir(fullfile(FC_dir,[namesInclSubj{iter_subj} '_419ROIcorrMat.mat']));
    foo = load(fullfile(FC_dir,RSfile.name)); 
    myCM(:,:) = foo.subj_final_corr_mat;
    X_nonreg(iter_subj,:) = jUpperTriMatToVec(myCM,1)'; 
end

nRois = size(myCM,1);
clear RSfile foo

%% 2. Regress out confounds from data

disp('(2) Regressing out confounds from data');

X_reg = CBIG_VK2019_regrOutConfounds(X_nonreg,confounds,behav_dir,motion_dir,RSfile_stem,commonInclSubj,namesInclSubj);
Y_reg = CBIG_VK2019_regrOutConfounds(Y_nonreg,confounds,behav_dir,motion_dir,RSfile_stem,commonInclSubj,namesInclSubj);

% Save original matrices
X0 = X_reg; Y0 = Y_reg;

%% 3. PLS analysis

disp('(3) Running PLS');

[U,S,V,Lx,Ly,explCovLC,LC_behav_loadings,LC_RSFC_loadings] = ...
    myPLS_analysis(X0,Y0,normalization_img,normalization_behav);

save(fullfile(out_dir,'PLSresults.mat'));

%% 4. Permutation testing

disp('(4) Permutation testing over LCs');

pvals_LC = myPLS_permut(X0,Y0,U,S,nPerms_rest,diagnosis_grouping,normalization_img,normalization_behav,1000);

% FDR correction over the first 5 LCs
[signif_LC, ~] = FDR(pvals_LC(1:5), 0.05);

% Display significant LCs
disp('Significant LCs (after FDR correction):');
for iter_lc = 1:length(signif_LC)
    this_lc = signif_LC(iter_lc);
    disp(['LC' num2str(this_lc) ' (p = ' num2str(pvals_LC(this_lc),'%0.3f') ') explains ' ...
        num2str(round(100*explCovLC(this_lc))) '% of covariance']);
end

save(fullfile(out_dir,['PLSresults_' num2str(nPerms_rest) 'permuts.mat']),...
    'pvals_LC','X0','Y0','U','S','nPerms_rest','diagnosis_grouping','normalization_img','normalization_behav');

%% 5. Plots

disp('(5) Plotting RSFC & behavior loadings and composite scores');

% RSFC & behavior loadings
CBIG_VK2019_plot_loadings(out_dir,LC_behav_loadings,behavNames,LC_RSFC_loadings,nRois,signif_LC); 

% RSFC & behavior composite scores
CBIG_VK2019_plot_subjScores(Lx,Ly,CONST_DIAGNOSIS,diagnosis_grouping,signif_LC); 

%% 6. Bootstrap RSFC & behavior loadings

disp('(6) Bootstrapping over RSFC & behavior loadings');

% Re-compute loadings using bootstrap 
[LC_RSFC_loadings_boot,LC_behav_loadings_boot,all_boot_orders] = CBIG_VK2019_bootstrap_loadings...
    (X0,Y0,U,signif_LC,nBootstraps,diagnosis_grouping,normalization_img,normalization_behav,1000);

save(fullfile(out_dir,['PLS_bootstrapLoadings_' num2str(nBootstraps) 'bootstraps.mat']),...
    'LC_RSFC_loadings_boot','LC_behav_loadings_boot','X0','Y0','U','signif_LC',...
    'nBootstraps','diagnosis_grouping','normalization_img','normalization_behav','all_boot_orders');

% Compute confidence intervals, z-scores, p-values of loadings
[std_behav_boot,zscore_behav_boot,pvals_behav_boot,std_RSFC_boot,zscore_RSFC_boot,pvals_RSFC_boot] = ...
    CBIG_VK2019_bootstrap_stats(LC_behav_loadings,LC_behav_loadings_boot,LC_RSFC_loadings,LC_RSFC_loadings_boot,...
    nRois,signif_LC,scripts_dir,out_dir);

save(fullfile(out_dir,['PLS_bootstrapResults_' num2str(nBootstraps) 'bootstraps.mat']),...
    'std_behav_boot','zscore_behav_boot','pvals_behav_boot','std_RSFC_boot','zscore_RSFC_boot','pvals_RSFC_boot',...
    'LC_behav_loadings','LC_behav_loadings_boot','LC_RSFC_loadings','LC_RSFC_loadings_boot','nRois',...
    'signif_LC','scripts_dir','out_dir');

%% 7. Cross-validation

disp('(7) Cross-validation');

nFolds = 5;
[corr_LxLy_tr,corr_LxLy_te,pvals_corr_LxLy_te] = CBIG_VK2019_crossval...
    (X0,Y0,signif_LC,diagnosis_grouping,nFolds,nPerms_cv);

save(fullfile(out_dir,['PLS_' num2str(nFolds) 'fold_crossval_' num2str(nPerms_cv) 'permuts.mat']),...
    'corr_LxLy_tr','corr_LxLy_te','pvals_corr_LxLy_te','X0','Y0','signif_LC','diagnosis_grouping','nFolds','nPerms_cv');

%% 8. Validation on task FC

disp('(8) Task FC validation');

% Project task FC onto RSFC saliences 
[corr_LxLy_task,pvals_corr_LxLy_task] = CBIG_VK2019_validation_taskFC...
    (Y_reg,U,V,taskFC_dir,signif_LC,nRois,namesInclSubj,...
    diagnosis_grouping,nPerms_task,normalization_img);

save(fullfile(out_dir,['PLS_taskValidation_' num2str(nPerms_task) 'permuts.mat']),...
    'corr_LxLy_task','pvals_corr_LxLy_task','Y_reg','U','V','taskFC_dir','signif_LC','nRois',...
    'namesInclSubj','commonInclSubj','diagnosis_grouping','nPerms_task','normalization_img');


rmpath(fullfile(getenv('CBIG_CODE_DIR'),'external_packages','matlab','non_default_packages','PLS_MIPlab'));
rmpath(scripts_dir);

