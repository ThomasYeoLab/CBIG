function decline = CBIG_hypoTest_3t_withinFactorAcrossStages(stats)

% decline = CBIG_hypoTest_3t_withinFactorAcrossStages(stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

beta = stats.Bhat;

betaCov = stats.CovBhat;

decline = zeros(9, 3);




fprintf('\n\n****************** Decline: Within Factor, Across Stages ******************\n\n');


%% ------------------------- Temporal

fprintf('\n----------------------------- Temporal -----------------------------\n\n');

% CN: beta_t = 0
C = [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('CN slope > 0, p = %e\n', F_C.pval);
else
    fprintf('CN slope < 0, p = %e\n', F_C.pval);
end
idx = find(C==1);
decline(1, :) = [beta(idx), sqrt(betaCov(idx, idx)), F_C.pval];

% Overall: beta_mt = beta_dt = 0
C = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall p = %e\n', F_C.pval);

% MCI - CN: beta_mt = 0
C = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('CN steeper than MCI, p = %e\n', F_C.pval);
else
    fprintf('MCI steeper than CN, p = %e\n', F_C.pval);
end
idx = find(C==1);
decline(2, :) = [beta(idx), sqrt(betaCov(idx, idx)), F_C.pval];

% AD - MCI: beta_dt - beta_mt = 0
C = [0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('MCI steeper than AD, p = %e\n', F_C.pval);
else
    fprintf('AD steeper than MCI, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
decline(3, :) = [...
    beta(idx1)-beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)),...
    F_C.pval];







%% ------------------------- Subcortical 

fprintf('\n----------------------------- Subcortical -----------------------------\n\n');

% CN: beta_t + beta_st = 0
C = [0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('CN slope > 0, p = %e\n', F_C.pval);
else
    fprintf('CN slope < 0, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(4, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% Overall: beta_mt + beta_mst = beta_dt + beta_dst = 0
C = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall p = %e\n', F_C.pval);

% MCI - CN: beta_mt + beta_mst = 0
C = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('CN steeper than MCI, p = %e\n', F_C.pval);
else
    fprintf('MCI steeper than CN, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(5, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% AD - MCI: beta_dt + beta_dst - beta_mt - beta_mst = 0
C = [0 0 0 0 0 0 0 0 0 0 -1 1 0 0 -1 0 1 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('MCI steeper than AD, p = %e\n', F_C.pval);
else
    fprintf('AD steeper than MCI, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(6, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];







%% ------------------------- Cortical

fprintf('\n----------------------------- Cortical -----------------------------\n\n');

% CN: beta_t + beta_ct = 0
C = [0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('CN slope > 0, p = %e\n', F_C.pval);
else
    fprintf('CN slope < 0, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(7, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% Overall: beta_mt + beta_mct = beta_dt + beta_dct = 0
C = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall p = %e\n', F_C.pval);

% MCI - CN: beta_mt + beta_mct = 0
C = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('CN steeper than MCI, p = %e\n', F_C.pval);
else
    fprintf('MCI steeper than CN, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(8, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% AD - MCI: beta_dt + beta_dct - beta_mt - beta_mct = 0
C = [0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 -1 0 1 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0
    fprintf('MCI steeper than AD, p = %e\n', F_C.pval);
else
    fprintf('AD steeper than MCI, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(9, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];
