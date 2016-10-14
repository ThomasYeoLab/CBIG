function [baseline, decline] = CBIG_hypoTest_3t(stats)

% [baseline, decline] = CBIG_hypoTest_3t(stats)
beta = stats.Bhat;
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

betaCov = stats.CovBhat;

baseline = zeros(9, 3);
decline = zeros(9, 3);






%% ------------------------- AD

fprintf('\n----------------------------- AD -----------------------------\n');

%% -------- Baseline

fprintf('\nBaseline\n\n');

% Overall: beta_s + beta_ds = beta_c + beta_dc = 0
C = [0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_s + beta_ds = 0
C = [0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
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

% C - T: beta_c + beta_dc = 0
C = [0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c + beta_dc > 0
    fprintf('T < C, p = %e\n', F_C.pval);
else
    fprintf('T > C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(2, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% C - S: beta_s + beta_ds = beta_c + beta_dc
C = [0 0 0 -1 1 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c + beta_dc > beta_s + beta_ds
    fprintf('S < C, p = %e\n', F_C.pval);
else
    fprintf('S > C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(3, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% Overall: beta_st + beta_dst = beta_ct + beta_dct = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_st + beta_dst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0];
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

% C - T: beta_ct + beta_dct = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct + beta_dct > 0
    fprintf('T > C, p = %e\n', F_C.pval);
else
    fprintf('T < C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(2, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% C - S: beta_ct + beta_dct = beta_st + beta_dst
C = [0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 -1 1 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct + beta_dct > beta_st + beta_dst
    fprintf('S > C, p = %e\n', F_C.pval);
else
    fprintf('S < C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(C==-1);
idx3 = ind(1);
idx4 = ind(2);
decline(3, :) = [...
    beta(idx1)+beta(idx2)-beta(idx3)-beta(idx4),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+betaCov(idx3, idx3)+betaCov(idx4, idx4)+2*betaCov(idx1, idx2)-2*betaCov(idx1, idx3)-2*betaCov(idx1, idx4)-2*betaCov(idx2, idx3)-2*betaCov(idx2, idx4)+2*betaCov(idx3, idx4)),...
    F_C.pval];






%% ------------------------- MCI

fprintf('\n----------------------------- MCI -----------------------------\n');

%% -------- Baseline

fprintf('\nBaseline\n\n');

% Overall: beta_s + beta_ms = beta_c + beta_mc = 0
C = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_s + beta_ms = 0
C = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_s + beta_ms > 0
    fprintf('T < S, p = %e\n', F_C.pval);
else
    fprintf('T > S, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(4, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% C - T: beta_c + beta_mc = 0
C = [0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c + beta_mc > 0
    fprintf('T < C, p = %e\n', F_C.pval);
else
    fprintf('T > C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(5, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% C - S: beta_s + beta_ms = beta_c + beta_mc
C = [0 0 0 -1 1 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c + beta_mc > beta_s + beta_ms
    fprintf('S < C, p = %e\n', F_C.pval);
else
    fprintf('S > C, p = %e\n', F_C.pval);
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

% Overall: beta_st + beta_mst = beta_ct + beta_mct = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_st + beta_mst = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_st + beta_mst > 0
    fprintf('T > S, p = %e\n', F_C.pval);
else
    fprintf('T < S, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(4, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% C - T: beta_ct + beta_mct = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct + beta_mct > 0
    fprintf('T > C, p = %e\n', F_C.pval);
else
    fprintf('T < C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(5, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

% C - S: beta_st + beta_mst = beta_ct + beta_mct
C = [0 0 0 0 0 0 0 0 0 0 0 0 -1 1 -1 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct + beta_mct > beta_st + beta_mst
    fprintf('S > C, p = %e\n', F_C.pval);
else
    fprintf('S < C, p = %e\n', F_C.pval);
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






%% ------------------------- Normal

fprintf('\n----------------------------- Normal -----------------------------\n');

%% Baseline

fprintf('\nBaseline\n\n');

% Overall: beta_s = beta_c = 0
C = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_s = 0
C = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_s > 0
    fprintf('T < S, p = %e\n', F_C.pval);
else
    fprintf('T > S, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
baseline(7, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% C - T: beta_c = 0
C = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c > 0
    fprintf('T < C, p = %e\n', F_C.pval);
else
    fprintf('T > C, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
baseline(8, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% C - S: beta_s = beta_c
C = [0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c > beta_s
    fprintf('S < C, p = %e\n', F_C.pval);
else
    fprintf('S > C, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
baseline(9, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% Overall: beta_st = beta_ct = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
fprintf('Overall: p = %e\n', F_C.pval);

% S - T: beta_st = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_st > 0
    fprintf('T > S, p = %e\n', F_C.pval);
else
    fprintf('T < S, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
decline(7, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% C - T: beta_ct = 0
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct > 0
    fprintf('T > C, p = %e\n', F_C.pval);
else
    fprintf('T < C, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
decline(8, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

% C - S: beta_st = beta_ct
C = [0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct > beta_st
    fprintf('S > C, p = %e\n', F_C.pval);
else
    fprintf('S < C, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
idx2 = find(C==-1);
decline(9, :) = [beta(idx1)-beta(idx2), sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)-2*betaCov(idx1, idx2)), F_C.pval];
