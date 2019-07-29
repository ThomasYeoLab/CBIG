function [corr_LxLy_tr,corr_LxLy_te,pvals_corr_LxLy_te] = CBIG_VK2019_crossval(X,Y,signif_LC,grouping,k,nPerms)

% This function computes PLS in k-fold cross-validation.
% The sample is first divided into training & test set, accounting for the
% diagnostic groups. 
% PLS is computed on the training data to obtain the saliences, which are 
% then projected on the test data.
% Correlations betwen RSFC & behavior composite scores for significant 
% latent components are computed for the training and test sets.
% Permutation testing is used to obtain p-values for
% the RSFC-behavior correlations for each test fold. 
%
% Inputs:
% - X                   : N x M matrix, N is #subjects, M is #FC, RSFC data 
% - Y                   : N x B matrix, N is #subjects, B is #behaviors, behavior data
% - signif_LC           : significant latent components (e.g. [1,2])
% - grouping            : N x 1 vector, subject (diagnostic) grouping 
%                         e.g. grouping = [1,2,3]
% - k                   : number of folds for cross-validation
% - nPerms              : number of permutations
%
% Outputs:
% - corr_LxLy_tr        : correlations between RSFC and behavior composite
%                         scores of training data
% - corr_LxLy_te        : correlations between RSFC and behavior composite
%                         scores of testing data
% - pvals_corr_LxLy_te  : p-values of correlations between RSFC and 
%                         behavior composite scores of testing data
%
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rng('default');

% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nSubj = size(X,1);
nGroups = size(unique(grouping),1);
subj_grouping = ones(nSubj,1);

% Set subjects within folds        
subj_nums = double(1:nSubj)';

for iter_group = 1:nGroups
    
    nSubj_groups(iter_group) = size(find(grouping==iter_group),1);
    
    C = cvpartition(nSubj_groups(iter_group),'Kfold',k);
    nSubj_groups_tr(iter_group,:) = C.TrainSize;
    nSubj_groups_te(iter_group,:) = C.TestSize;
    
    these_subj{iter_group} = subj_nums(find(grouping==iter_group));
    
    a=1;
    for iter_fold = 1:k
        these_subj_te{iter_group,iter_fold} = these_subj{iter_group}(a : a -1 + nSubj_groups_te(iter_group,iter_fold));
        these_subj_tr{iter_group,iter_fold} = setdiff(these_subj{iter_group},these_subj_te{iter_group,iter_fold});
        a = a + nSubj_groups_te(iter_group,iter_fold);
    end
    clear a C
end
        

%%% 1. PLS on training sample, then project saliences of training sets on 
% test sets 

for iter_fold = 1:k    
           
    idx_tr = []; idx_te = [];
    
    for iter_group = 1:nGroups
        idx_tr = [idx_tr ; these_subj_tr{iter_group,iter_fold}];
        idx_te = [idx_te ; these_subj_te{iter_group,iter_fold}];
    end
        
    % Separate X and Y into training and testing   
    X_tr = X(idx_tr,:);
    X_te = X(idx_te,:);
    
    Y_tr = Y(idx_tr,:);
    Y_te = Y(idx_te,:);
    
    subj_grouping_tr = subj_grouping(idx_tr);   
    all_grouping_te{iter_fold} = grouping(idx_te);
    
    % Normalize X_tr and Y_tr    
    X_tr = X_tr;
    Y_tr = Y_tr;
    
    % Z-score manually to save mean & std for later
    [X_tr, MeanX, SigmaX] = zscore(X_tr);
    [Y_tr, MeanY, SigmaY] = zscore(Y_tr);
    
    allX_tr{iter_fold} = X_tr;
    allY_tr{iter_fold} = Y_tr;
    
    % Cross-covariance matrix
    R_tr = myPLS_cov(X_tr,Y_tr,1,subj_grouping_tr);
    
    % Singular value decomposition   
    [U_tr, S_tr, V_tr] = svd(R_tr,'econ');
    nLCs = min(size(S_tr));
    
    % ICA convention: turn LCs such that max is positive
    for iter_lc = 1:nLCs
        [~,idx] = max(abs(V_tr(:,iter_lc)));
        if sign(V_tr(idx,iter_lc))<0
            V_tr(:,iter_lc) = -V_tr(:,iter_lc);
            U_tr(:,iter_lc) = -U_tr(:,iter_lc);
        end
    end
    
    allU_tr{iter_fold} = U_tr ;
    allV_tr{iter_fold} = V_tr ;
    
    % RSFC & behavior composite scores (training data)    
    Lx_tr = X_tr * V_tr;
    Ly_tr = Y_tr * U_tr;
    
    allLx_tr{iter_fold} = Lx_tr;
    allLy_tr{iter_fold} = Ly_tr;
    
    % Normalize X_te and Y_te using mean & sigma from training data    
    sigma0 = SigmaX;
    sigma0(sigma0==0) = 1;
    X_te = bsxfun(@minus,X_te, MeanX);
    X_te = bsxfun(@rdivide, X_te, sigma0);
    
    sigma0 = SigmaY;
    sigma0(sigma0==0) = 1;
    Y_te = bsxfun(@minus,Y_te, MeanY);
    Y_te = bsxfun(@rdivide, Y_te, sigma0);
    
    allX_te{iter_fold} = X_te;
    allY_te{iter_fold} = Y_te;
       
    % RSFC & behavior composite scores (testing data)
    Lx_te = X_te * V_tr;
    Ly_te = Y_te * U_tr;
    
    allLx_te{iter_fold} = Lx_te;
    allLy_te{iter_fold} = Ly_te;
    
    % Correlations between RSFC and behavior composite scores
    % (training & testing data) within each fold
    for iter_lc = 1:size(signif_LC,1)
        corr_LxLy_tr(iter_fold,iter_lc) = corr(allLx_tr{iter_fold}(:,iter_lc),allLy_tr{iter_fold}(:,iter_lc));
        corr_LxLy_te(iter_fold,iter_lc) = corr(allLx_te{iter_fold}(:,iter_lc),allLy_te{iter_fold}(:,iter_lc));
    end
    
    %%% 2. Permutation testing on RSFC-behavior correlations across folds

    for iter_perm = 1:nPerms
        
        Yp_te = [];
        for iter_group = 1:nGroups
            thisY_te = Y_te(find(all_grouping_te{iter_fold} == iter_group),:);
            perm_order = randperm(size(thisY_te,1));
            thisYp_te = thisY_te(perm_order,:);
            Yp_te = [Yp_te; thisYp_te];
            clear thisY_te thisYp_te
        end
        
        Lyp_te = Yp_te * U_tr;
        allLy_p_te{iter_fold} = Lyp_te;
        
        for iter_lc = 1:size(signif_LC,1)
            corr_LxLyp_te(iter_fold,iter_lc) = corr(allLx_te{iter_fold}(:,iter_lc),allLy_p_te{iter_fold}(:,iter_lc));
        end
        
    end
    
    % Compute p-value for correlations in test sets
    for iter_lc = 1:size(signif_LC,1)
        pvals_corr_LxLy_te(iter_fold,iter_lc) = (sum(corr_LxLyp_te(:,iter_lc) > ...
            corr_LxLy_te(iter_fold,iter_lc)) +1)/(nPerms+1) ;
    end 

end

% Display mean correlations between RSFC & behavior composite scores
 disp('Correlations between RSFC & behavior composite scores for training & testing data:');

for iter_lc = 1:size(signif_LC,1)

    disp(['LC' num2str(iter_lc)]);

    this_corr = corr_LxLy_tr(:,iter_lc);
    disp(['Training data: Mean (SD) corr Lx-Ly = ' num2str(mean(this_corr),'%0.2f') ...
        ' (' num2str(std(this_corr),'%0.2f') ') ' num2str(min(this_corr),'%0.2f') '-' ...
        num2str(max(this_corr),'%0.2f')]);
    clear this_corr

    this_corr = corr_LxLy_te(:,iter_lc);
    disp(['Test data: Mean (SD) corr Lx-Ly = ' num2str(mean(this_corr),'%0.2f') ...
        ' (' num2str(std(this_corr),'%0.2f') ') ' num2str(min(this_corr),'%0.2f') '-' ...
        num2str(max(this_corr),'%0.2f')]);
    clear this_corr

    disp(['Range permuted p = ' num2str(min(pvals_corr_LxLy_te(:,iter_lc)),'%0.3f') ' - ' ...
        num2str(max(pvals_corr_LxLy_te(:,iter_lc)),'%0.3f')]);

end

