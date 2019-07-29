function [corr_LxLy_task,pvals_corr_LxLy_task] = CBIG_VK2019_validation_taskFC...
    (Y_rest,U_rest,V_rest,taskFC_dir,signif_LC,nRois,names_subj,grouping,nPerms,normalization_img)
% 
% This function aims to validate PLS latent components (LCs) obtained with
% RSFC on task FC. The RSFC saliences are projected onto the task FC data 
% to obtain task FC composite scores. Correlations between task FC & 
% behavior composite scores (which are the same as for the RSFC model) are
% then computed for all significant LCs. 
% Permutation testing (accounting for diagnostic groups) is used to 
% compute pvalues for the correlations between task FC & behavior composite 
% scores.
%
% Inputs:
% - Y_rest                   : N x B matrix, N is #subjects, B is
%                              #behaviors, RSFC data
% - U_rest                   : B x L matrix, L is #LCs, behavior saliences
%                              of RSFC PLS model
% - V_rest                   : M x L matrix, M is # FC, RSFC saliences of
%                              RSFC PLS model
% - taskFC_dir               : path to directory where task FC data 
%                              are located (as correlation matrices saved as
%                              '<name_subject>_<nRois>_ROIcorrMat.mat')
% - signif_LC                : significant latent components (e.g. [1,2])
% - nRois                    : number of ROIs in RSFC matrix (e.g. 419)
% - names_subj               : N x 1 string, names of subjects
% - grouping                 : N x 1 vector, subject (diagnostic) grouping 
%                              e.g. grouping = [1,2,3]
% - nPerms                   : number of permutations
% - normalization_img        : normalization options for task FC data
% - normalization_behav      : normalization options for behavior data
%                              0 = no normalization
%                              1 = zscore across all subjects
%                              2 = zscore within groups (default)
%                              3 = std normalization across subjects (no centering)
%                              4 = std normalization within groups (no centering)
%
% Outputs:
% - corr_LxLy_task           : correlation between task FC and behavior 
%                              composite scores
% - pvals_corr_Lx_Ly_task    : p-values of correlations between task FC and
%                              behavior composite scores
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

nSubj = numel(names_subj);
nGroups = size(unique(grouping),1);
nFC = size(V_rest,1);

% Load task FC
taskFC_vec = zeros(nSubj,nFC);
missing = []; 

for iter_subj = 1:nSubj    
    thisFile = fullfile(taskFC_dir,[names_subj{iter_subj} '_' num2str(nRois) 'ROIcorrMat.mat']);
    if exist(thisFile)
        load(thisFile);
        taskFC_vec(iter_subj,:) = jUpperTriMatToVec(subj_final_corr_mat,1); clear subj_final_corr_mat
    else
        taskFC_vec(iter_subj,:) = zeros(1,size(taskFC_vec,2));
        disp([names_subj{iter_subj} ' - no task FC']);
        missing = [missing iter_subj];
    end
    
end

disp([num2str(size(missing,2)) ' subjects missing - excluded from task FC validation']);

taskFC_vec(isnan(taskFC_vec)) = 0;

% Remove subjects with missing data
nSubj = nSubj - size(missing,2);
taskFC_vec(missing,:) = [];
Y_rest(missing,:) = [];
names_subj(missing) = [];
grouping(missing) = [];

% Check that dimensions of taskFC_vec & Y are correct
if(size(taskFC_vec,1) ~= size(Y_rest,1))
    error('Input arguments taskFC_vec and Y_rest should have the same number of rows');
end

% Normalize task FC (Y is already normalized)
subj_grouping = ones(nSubj,1);
X_task = myPLS_norm(taskFC_vec,1,subj_grouping,normalization_img);

% Compute task FC composite scores 
Lx_task = X_task * V_rest;
Ly_task = Y_rest * U_rest; % Ly_rest;

% Compute correlations between task FC & behavior composite scores
for iter_lc = 1:size(signif_LC,1)
    corr_LxLy_task(iter_lc) = corr(Lx_task(:,iter_lc),Ly_task(:,iter_lc));
end

%% Permutation testing 
% over correlations between RSFC & behavioral composite scores
       
disp('... Permutations ...')
for iter_perm = 1:nPerms
    
    % Display number of permutations (every 50 permuts)
    if mod(iter_perm,50)==0, disp(num2str(iter_perm)); end
    
    % X is already normalized
    Xp = X_task;
    
    % Permute (normalized) Y within each diagnostic group
    Yp = []; % Lyp
    for iter_group = 1:nGroups
        clear thisY0 thisYp
        thisY = Y_rest(find(grouping==iter_group),:);
        perm_order = randperm(size(thisY,1));
        thisYp = thisY(perm_order,:);
        Yp = [Yp; thisYp];
    end
            
    % Behavior composite scores under null hypothesis
    Lyp_task = Yp * U_rest;   
       
    % Compute correlation between task FC composite scores & behavior 
    % composite scores obtained by permuting rows of Y_rest
    for iter_lc = 1:size(signif_LC,1)
        [corr_LxLyp_task(iter_perm,iter_lc),~] = corr(Lx_task(:,iter_lc),Lyp_task(:,iter_lc));              
    end
    
end

% Compute & display significance of correlations
disp('Correlations between task FC and behavioral composite scores');
for iter_lc = 1:size(signif_LC,1)
    a = size(find(corr_LxLyp_task(:,iter_lc) > corr_LxLy_task(iter_lc)),1);
    pvals_corr_LxLy_task(iter_lc,1) = (a + 1)/(nPerms + 1) ; clear a
    disp(['LC' num2str(iter_lc) ' r = ' num2str(corr_LxLy_task(iter_lc),'%0.2f') ....
        ' permuted p = ' num2str(pvals_corr_LxLy_task(iter_lc),'%0.3f')]);
end

