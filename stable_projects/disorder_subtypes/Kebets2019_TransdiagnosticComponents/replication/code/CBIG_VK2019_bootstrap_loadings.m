
function [LC_RSFC_loadings_boot,LC_behav_loadings_boot,all_boot_orders] = CBIG_VK2019_bootstrap_loadings...
    (X,Y,U,signif_LC,nBootstraps,grouping,normalization_img,normalization_behav,seed)
% 
% This function computes bootstrap resampling with replacement on X and Y, 
% accounting for diagnostic groups.
% The PLS analysis is re-computed for each bootstrap sample. 
% The RSFC & behavior loadings are also computed for each bootstrap sample
% for all significant latent components. 
%
% Inputs:
% - X                        : N x M matrix, N is #subjects, M is #FC, RSFC data 
% - Y                        : N x B matrix, B is #behaviors, behavior data
% - U                        : B x L matrix, L is # latent components (LCs),
%                              behavior saliences
% - signif_LC                : significant LCs
% - nBootstraps              : number of bootstrap samples
% - grouping                 : N x 1 vector, subject (diagnostic) grouping 
% - normalization_img        : normalization options for FC data
% - normalization_behav      : normalization options for behavior data
%                              0 = no normalization
%                              1 = zscore across all subjects
%                              2 = zscore within groups (default)
%                              3 = std normalization across subjects (no centering)
%                              4 = std normalization within groups (no centering)
% - seed                     : seed for random number generator, e.g. 1
%                             (default is 1000)
% 
% Outputs:
% - LC_behav_loadings_boot   : B x S x P matrix, S is # significant LCs, P is #bootstrap samples,
%                              bootstrapped behavior loadings for significant LCs 
% - LC_RSFC_loadings_boot    : M x S x P matrix, bootstrapped RSFC loadings
%                              for significant LCs
% - all_boot_orders          : N x P, subject index in all bootstrap samples
%
% Written by Valeria Kebets & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% Set random number generator 
if isempty(seed)
    rng(1000);
end

% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nSubj = size(X,1);
groupIDs = unique(grouping);
nGroups = length(groupIDs);
subj_grouping = ones(nSubj,1);

% Get bootstrap subject sampling
all_boot_orders = nan(nSubj,nBootstraps);

for iter_group = 1:nGroups
    
   % Get IDs of subjects in this group
   groupID_idx = find(grouping == groupIDs(iter_group));
   nSubj_group = length(groupID_idx);

   % Compute bootstrapping orders for the number of subjects in this group
   [boot_order_tmp,~] = rri_boot_order(nSubj_group,1,nBootstraps);

   % Change the order of the group indices according to the bootstrapping
   all_boot_orders(groupID_idx,:) = groupID_idx(boot_order_tmp);
end
        
        
for iter_boot = 1:nBootstraps
    
    if mod(iter_boot,50)==0, disp(num2str(iter_boot)); end
    
    % Resampling X
    Xb = X(all_boot_orders(:,iter_boot),:);
    Xb = myPLS_norm(Xb,1,subj_grouping,normalization_img);
    
    % Resampling Y
    Yb = Y(all_boot_orders(:,iter_boot),:);
    Yb = myPLS_norm(Yb,1,subj_grouping,normalization_behav);
    
    % Cross-covariance matrix between resampled X and Y
    Rb = myPLS_cov(Xb,Yb,1,subj_grouping);
    
    % Singular value decomposition of Rp
    [Ub,~,Vb] = svd(Rb,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Ub);
    Vb = Vb * rotatemat;
    Ub = Ub * rotatemat;
    
    % RSFC & behavior composite scores
    Lxb = Xb * Vb;
    Lyb = Yb * Ub;
    
    % Loadings   
    for iter_lc = 1:size(signif_LC,1)
        this_lc = signif_LC(iter_lc);
        
        % RSFC
        for iter_fc = 1:size(Xb,2)
            tmpy = Lxb(:,this_lc);
            tmpx = Xb(:,iter_fc);
            [r,~] = corrcoef(tmpx,tmpy.');
            these_RSFC_loadings(iter_fc,iter_lc) = r(1,2);
            clear tmpy tmpx r
        end
    
        % Behavior
        for iter_behav = 1:size(Y,2)
            tmpy = Lyb(:,this_lc);
            tmpx = Yb(:,iter_behav);
            [r,~] = corrcoef(tmpx,tmpy.');
            these_behav_loadings(iter_behav,iter_lc) = r(1,2);
            clear tmpy tmpx r
        end
    end
    
    LC_RSFC_loadings_boot(:,:,iter_boot) = these_RSFC_loadings;
    LC_behav_loadings_boot(:,:,iter_boot) = these_behav_loadings;      
    
end
     