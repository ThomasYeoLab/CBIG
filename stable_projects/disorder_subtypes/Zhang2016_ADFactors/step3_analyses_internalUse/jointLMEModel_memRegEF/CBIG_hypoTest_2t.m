function [baseline, decline] = CBIG_hypoTest_2t(stats)

% [baseline, decline] = CBIG_hypoTest_2t(stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

beta = stats.Bhat;

betaCov = stats.CovBhat;

baseline = zeros(3, 3);
decline = zeros(3, 3);






%% ------------------------- AD

fprintf('\n----------------------------- AD -----------------------------\n');

%% -------- Baseline

fprintf('\nBaseline\n\n');

% C - T+S: beta_c + beta_dc = 0
C = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c + beta_dc > 0
    fprintf('T+S < C, p = %e\n', F_C.pval);
else
    fprintf('T+S > C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(1, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% C - T+S: beta_ct + beta_dct = 0
C = [0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct + beta_dct > 0
    fprintf('T+S > C, p = %e\n', F_C.pval);
else
    fprintf('T+S < C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(1, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];






%% ------------------------- MCI

fprintf('\n----------------------------- MCI -----------------------------\n');

%% -------- Baseline

fprintf('\nBaseline\n\n');

% C - T+S: beta_c + beta_mc = 0
C = [0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c + beta_mc > 0
    fprintf('T+S < C, p = %e\n', F_C.pval);
else
    fprintf('T+S > C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(2, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% C - T+S: beta_ct + beta_mct = 0
C = [0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct + beta_mct > 0
    fprintf('T+S > C, p = %e\n', F_C.pval);
else
    fprintf('T+S < C, p = %e\n', F_C.pval);
end
ind = find(C==1);
idx1 = ind(1);
idx2 = ind(2);
decline(2, :) = [...
    beta(idx1)+beta(idx2),...
    sqrt(betaCov(idx1, idx1)+betaCov(idx2, idx2)+2*betaCov(idx1, idx2)),...
    F_C.pval];






%% ------------------------- Normal

fprintf('\n----------------------------- Normal -----------------------------\n');

%% Baseline

fprintf('\nBaseline\n\n');

% C - T+S: beta_c = 0
C = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_c > 0
    fprintf('T+S < C, p = %e\n', F_C.pval);
else
    fprintf('T+S > C, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
baseline(3, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];

%% -------- Decline

fprintf('\nDecline rate\n\n');

% C - T+S: beta_ct = 0
C = [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
F_C = lme_F(stats, C);
if F_C.sgn > 0 % beta_ct > 0
    fprintf('T+S > C, p = %e\n', F_C.pval);
else
    fprintf('T+S < C, p = %e\n', F_C.pval);
end
idx1 = find(C==1);
decline(3, :) = [beta(idx1), sqrt(betaCov(idx1, idx1)), F_C.pval];
