function [U,S,V,Lx,Ly,explCovLC,LC_behav_loadings,LC_RSFC_loadings] = myPLS_analysis(X,Y,normalization_img,normalization_behav)
% PLS analysis (main script)
%
% Inputs:
% - X                   : N x M matrix, N is #subjects, M is #FC 
% - Y                   : N x B matrix, N is #subjects, B is #behaviors
% - normalization_img   : normalization options for FC data
% - normalization_behav : normalization options for behavior data
%                         0 = no normalization
%                         1 = zscore across all subjects
%                         2 = zscore within groups (default)
%                         3 = std normalization across subjects (no centering)
%                         4 = std normalization within groups (no centering)
%
% Outputs:
% - U                   : B x L matrix, behavior saliences
% - S                   : L x L matrix, singular values
% - V                   : M x L matrix, RSFC saliences
% - Lx                  : N x L matrix, RSFC composite scores
% - Ly                  : N x L matrix, behavior composite scores
% - explCovLC           : covariance explained by each latent component
% - LC_behav_loadings   : B x L matrix, behavior loadings
% - LC_RSFC_loadings    : M x L matrix, RSFC loadings
%
% Written by Valeria Kebets and the MIPlab, with subfunctions borrowed
% from PLS toolbox by Rotman Baycrest
% (https://www.rotman-baycrest.on.ca/index.php?section=84)
%
% Please cite the following papers when using this code:
%
% Zoller D, Schaer M, Scariati E, Padula MC, Eliez S, Van De Ville D (2017).
% Disentangling resting-state BOLD variability and PCC
% functional connectivity in 22q11.2 deletion syndrom.
% Neuroimage 149, pp. 85-97.
%
% McIntosh AR, Lobaugh NJ (2004). Partial least squares analysis of
% neuroimaging data: applications and advances.
% Neuroimage 23(Suppl 1), pp. S250-263.


% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nSubj = size(X,1);
subj_grouping = ones(nSubj,1);

% Data normalization
X = myPLS_norm(X,1,subj_grouping,normalization_img);
Y = myPLS_norm(Y,1,subj_grouping,normalization_behav);

% Cross-covariance matrix
R = myPLS_cov(X,Y,1,subj_grouping);

% Singular value decomposition
[U,S,V] = svd(R,'econ');
nLCs = min(size(S)); % Number of latent components (LC)

% ICA convention: turn LCs such that max is positive
for iter_lc = 1:nLCs
    [~,idx] = max(abs(V(:,iter_lc)));
    if sign(V(idx,iter_lc))<0
        V(:,iter_lc) = -V(:,iter_lc);
        U(:,iter_lc) = -U(:,iter_lc);
    end
end

% Amount of covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));

% RSFC & behavioral composite scores
Lx = X * V;
Ly = Y * U;

% RSFC loadings (Pearson's correlations between Lx and X) 
for iter_lc = 1:nLCs
    for iter_fc = 1:size(X,2)
        clear tmpy tmpx r
        tmpy = Lx(:,iter_lc);
        tmpx = X(:,iter_fc);
        r = corrcoef(tmpx,tmpy.');
        LC_RSFC_loadings(iter_fc,iter_lc) = r(1,2);
    end
end

% Behavioral loadings (Pearson's correlations between Ly and Y) 
for iter_lc = 1:nLCs
    for iter_behav = 1:size(Y,2)
        clear tmpy tmpx r
        tmpy = Ly(:,iter_lc);
        tmpx = Y(:,iter_behav);
        r = corrcoef(tmpx,tmpy.');
        LC_behav_loadings(iter_behav,iter_lc) = r(1,2);
    end
end