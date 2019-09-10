function pvals_LC = myPLS_permut(X,Y,U,S,nPerms,grouping,normalization_img,normalization_behav,seed)
% Permutation testing over singular values obtained with PLS
% Rows (subjects) of Y are permuted within each diagnostic group
%
% Inputs:
% - X                    : N x M matrix, N is #subjects, M is #FC 
% - Y                    : N x B matrix, N is #subjects, B is #behaviors
% - U                    : B x L matrix, behavior saliences
% - S                    : L x L matrix, singular values
% - nPerms               : number of permutations
% - grouping             : N x 1 vector, subject grouping 
% - normalization_img    : normalization options for FC data
% - normalization_behav  : normalization options for behavior data
%                          0 = no normalization
%                          1 = zscore across all subjects
%                          2 = zscore within groups (default)
%                          3 = std normalization across subjects (no centering)
%                          4 = std normalization within groups (no centering)
% - seed                 : seed for random number generator, e.g. 1
%                          (default is 1000)
%
% Outputs:
% - pvals_LC             : p-values for each latent component
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


% Set random number generator 
if isempty(seed)
    rng(1000);
end


% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nSubj = size(X,1);
nGroups = size(unique(grouping),1);
subj_grouping = ones(nSubj,1);

disp('... Permutations ...')
for iter_perm = 1:nPerms
    
    % Display number of permutations (every 50 permuts)
    if mod(iter_perm,50)==0, disp(num2str(iter_perm)); end
    
    % Normalization of X
    Xp = myPLS_norm(X,1,subj_grouping,normalization_img);
    
    % Permute Y within each diagnostic group
    Yp = [];
    for iter_group = 1:nGroups
        clear thisY thisYp
        thisY = Y(find(grouping==iter_group),:);
        perm_order = randperm(size(thisY,1));
        thisYp = thisY(perm_order,:);
        Yp = [Yp; thisYp];
    end
    
    % Normalization of Y
    Yp = myPLS_norm(Yp,1,subj_grouping,normalization_img);
    
    % Cross-covariance matrix between X and permuted Y
    Rp = myPLS_cov(Xp,Yp,1,subj_grouping);
    
    % Singular value decomposition of Rp
    [Up,Sp,Vp] = svd(Rp,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Up);
    Up = Up * Sp * rotatemat;
    Sp = sqrt(sum(Up.^2));
    
    permsamp(:,iter_perm) = Sp';
    
    if iter_perm == 1
        sp= (Sp' >= diag(S));
    else
        sp = sp + (Sp'>= diag(S));
    end
    
end

% Compute p-values for each LC
pvals_LC = (sp + 1) ./ (nPerms + 1);
