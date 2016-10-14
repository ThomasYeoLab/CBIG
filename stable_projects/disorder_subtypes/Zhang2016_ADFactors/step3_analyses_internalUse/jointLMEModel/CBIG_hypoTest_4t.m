function [baseline, decline] = CBIG_hypoTest_4t(stats)

% [baseline, decline] = CBIG_hypoTest_4t(stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

beta = stats.Bhat;

betaCov = stats.CovBhat;

baseline = zeros(18, 3);
decline = zeros(18, 3);






%% ------------------------- AD

fprintf('\n----------------------------- AD -----------------------------\n');

%% -------- Baseline

fprintf('\nBaseline\n\n');

% Overall: beta_s + beta_ds = beta_f + beta_df = beta_p + beta_dp = 0
C = [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_s + beta_ds = 0
C = [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_s + beta_ds > 0
    fprintf('T < S, p = %e\n', F_C.pval);
else
    fprintf('T > S, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(1, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - T: beta_f + beta_df = 0
C = [0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_f + beta_df > 0
    fprintf('T < F, p = %e\n', F_C.pval);
else
    fprintf('T > F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(2, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% P - T: beta_p + beta_dp = 0
C = [0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p + beta_dp > 0
    fprintf('T < P, p = %e\n', F_C.pval);
else
    fprintf('T > P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(3, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - S:  beta_f + beta_df - beta_s - beta_ds = 0
C = [0 0 0 -1 1 0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_f + beta_df - beta_s - beta_ds > 0
    fprintf('S < F, p = %e\n', F_C.pval);
else
    fprintf('S > F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(4, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - S: beta_p + beta_dp - beta_s - beta_ds = 0
C = [0 0 0 -1 0 1 0 0 0 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p + beta_dp - beta_s - beta_ds > 0
    fprintf('S < P, p = %e\n', F_C.pval);
else
    fprintf('S > P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(5, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - F: beta_p + beta_dp - beta_f - beta_df = 0
C = [0 0 0 0 -1 1 0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p + beta_dp - beta_f - beta_df > 0
    fprintf('F < P, p = %e\n', F_C.pval);
else
    fprintf('F > P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(6, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% Overall: beta_st + beta_dst = beta_ft + beta_dft = beta_pt + beta_dpt = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_st + beta_dst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_st + beta_dst > 0
    fprintf('T > S, p = %e\n', F_C.pval);
else
    fprintf('T < S, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(1, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - T: beta_ft + beta_dft = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ft + beta_dft > 0
    fprintf('T > F, p = %e\n', F_C.pval);
else
    fprintf('T < F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(2, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% P - T: beta_pt + beta_dpt = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt + beta_dpt > 0
    fprintf('T > P, p = %e\n', F_C.pval);
else
    fprintf('T < P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(3, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - S:  beta_ft + beta_dft - beta_st - beta_dst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 -1 1 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ft + beta_dft - beta_st - beta_dst > 0
    fprintf('S > F, p = %e\n', F_C.pval);
else
    fprintf('S < F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(4, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - S: beta_pt + beta_dpt - beta_st - beta_dst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0 0 -1 0 1 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt + beta_dpt - beta_st - beta_dst > 0
    fprintf('S > P, p = %e\n', F_C.pval);
else
    fprintf('S < P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(5, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - F: beta_pt + beta_dpt - beta_ft - beta_dft = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 -1 1 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt + beta_dpt - beta_ft - beta_dft > 0
    fprintf('F > P, p = %e\n', F_C.pval);
else
    fprintf('F < P, p = %e\n', F_C.pval);
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






%% ------------------------- MCI

fprintf('\n----------------------------- MCI -----------------------------\n');

%% -------- Baseline

fprintf('\nBaseline\n\n');

% Overall: beta_s + beta_ms = beta_f + beta_mf = beta_p + beta_mp = 0
C = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_s + beta_ms = 0
C = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_s + beta_ms > 0
    fprintf('T < S, p = %e\n', F_C.pval);
else
    fprintf('T > S, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(7, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - T: beta_f + beta_mf = 0
C = [0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_f + beta_mf > 0
    fprintf('T < F, p = %e\n', F_C.pval);
else
    fprintf('T > F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(8, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% P - T: beta_p + beta_mp = 0
C = [0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p + beta_mp > 0
    fprintf('T < P, p = %e\n', F_C.pval);
else
    fprintf('T > P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(9, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - S:  beta_f + beta_mf - beta_s - beta_ms = 0
C = [0 0 0 -1 1 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_f + beta_mf - beta_s - beta_ms > 0
    fprintf('S < F, p = %e\n', F_C.pval);
else
    fprintf('S > F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(10, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - S: beta_p + beta_mp - beta_s - beta_ms = 0
C = [0 0 0 -1 0 1 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p + beta_mp - beta_s - beta_ms > 0
    fprintf('S < P, p = %e\n', F_C.pval);
else
    fprintf('S > P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(11, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - F: beta_p + beta_mp - beta_f - beta_mf = 0
C = [0 0 0 0 -1 1 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p + beta_mp - beta_f - beta_mf > 0
    fprintf('F < P, p = %e\n', F_C.pval);
else
    fprintf('F > P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(12, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% Overall: beta_st + beta_mst = beta_ft + beta_mft = beta_pt + beta_mpt = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_st + beta_mst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_st + beta_mst > 0
    fprintf('T > S, p = %e\n', F_C.pval);
else
    fprintf('T < S, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(7, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - T: beta_ft + beta_mft = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ft + beta_mft > 0
    fprintf('T > F, p = %e\n', F_C.pval);
else
    fprintf('T < F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(8, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% P - T: beta_pt + beta_mpt = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt + beta_mpt > 0
    fprintf('T > P, p = %e\n', F_C.pval);
else
    fprintf('T < P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(9, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% F - S:  beta_ft + beta_mft - beta_st - beta_mst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 -1 1 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ft + beta_mft - beta_st - beta_mst > 0
    fprintf('S > F, p = %e\n', F_C.pval);
else
    fprintf('S < F, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(10, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - S: beta_pt + beta_mpt - beta_st - beta_mst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 1 -1 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt + beta_mpt - beta_st - beta_mst > 0
    fprintf('S > P, p = %e\n', F_C.pval);
else
    fprintf('S < P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(11, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

% P - F: beta_pt + beta_mpt - beta_ft - beta_mft = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 -1 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt + beta_mpt - beta_ft - beta_mft > 0
    fprintf('F > P, p = %e\n', F_C.pval);
else
    fprintf('F < P, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(12, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];






%% ------------------------- Normal

fprintf('\n----------------------------- Normal -----------------------------\n');

%% Baseline

fprintf('\nBaseline\n\n');

% Overall: beta_s = beta_f = beta_p = 0
C = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_s = 0
C = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_s > 0
    fprintf('T < S, p = %e\n', F_C.pval);
else
    fprintf('T > S, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
baseline(13, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% F - T: beta_f = 0
C = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_f > 0
    fprintf('T < F, p = %e\n', F_C.pval);
else
    fprintf('T > F, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
baseline(14, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% P - T: beta_p = 0
C = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p > 0
    fprintf('T < P, p = %e\n', F_C.pval);
else
    fprintf('T > P, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
baseline(15, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% F - S: beta_f - beta_s = 0
C = [0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_f - beta_s > 0
    fprintf('S < F, p = %e\n', F_C.pval);
else
    fprintf('S > F, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
baseline(16, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];

% P - S: beta_p - beta_s = 0
C = [0 0 0 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p - beta_s > 0
    fprintf('S < P, p = %e\n', F_C.pval);
else
    fprintf('S > P, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
baseline(17, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];

% P - F: beta_p - beta_f = 0
C = [0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_p - beta_f > 0
    fprintf('F < P, p = %e\n', F_C.pval);
else
    fprintf('F > P, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
baseline(18, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% Overall: beta_st = beta_ft = beta_pt = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_st = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_st > 0
    fprintf('T > S, p = %e\n', F_C.pval);
else
    fprintf('T < S, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
decline(13, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% F - T: beta_ft = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ft > 0
    fprintf('T > F, p = %e\n', F_C.pval);
else
    fprintf('T < F, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
decline(14, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% P - T: beta_pt = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt > 0
    fprintf('T > P, p = %e\n', F_C.pval);
else
    fprintf('T < P, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
decline(15, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% F - S: beta_ft - beta_st = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ft - beta_st > 0
    fprintf('S > F, p = %e\n', F_C.pval);
else
    fprintf('S < F, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
decline(16, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];

% P - S: beta_pt - beta_st = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt - beta_st > 0
    fprintf('S > P, p = %e\n', F_C.pval);
else
    fprintf('S < P, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
decline(17, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];

% P - F: beta_pt - beta_ft = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_pt - beta_ft > 0
    fprintf('F > P, p = %e\n', F_C.pval);
else
    fprintf('F < P, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
decline(18, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];
