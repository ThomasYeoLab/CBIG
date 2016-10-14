function baseline = CBIG_hypoTest_4t(betas, stats)

% baseline = CBIG_hypoTest_4t(betas, stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

baseline = zeros(18, 3);


%% AD

fprintf('\n----------------------------- AD -----------------------------\n');

fprintf('\nBaseline\n\n');

% Overall: beta_s + beta_ds = beta_f + beta_df = beta_p + beta_dp = 0
H = [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0];
c = [0; 0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall: p = %e\n', p);

% S - T: beta_s + beta_ds = 0
H = [0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_s + beta_ds > 0
    fprintf('T < S, p = %e\n', p);
else
    fprintf('T > S, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(1, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% F - T: beta_f + beta_df = 0
H = [0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_f + beta_df > 0
    fprintf('T < F, p = %e\n', p);
else
    fprintf('T > F, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(2, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% P - T: beta_p + beta_dp = 0
H = [0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p + beta_dp > 0
    fprintf('T < P, p = %e\n', p);
else
    fprintf('T > P, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(3, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% F - S:  beta_f + beta_df - beta_s - beta_ds = 0
H = [0 0 0 -1 1 0 0 0 0 -1 1 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_f + beta_df - beta_s - beta_ds > 0
    fprintf('S < F, p = %e\n', p);
else
    fprintf('S > F, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(4, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];

% P - S: beta_p + beta_dp - beta_s - beta_ds = 0
H = [0 0 0 -1 0 1 0 0 0 -1 0 1 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p + beta_dp - beta_s - beta_ds > 0
    fprintf('S < P, p = %e\n', p);
else
    fprintf('S > P, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(5, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];

% P - F: beta_p + beta_dp - beta_f - beta_df = 0
H = [0 0 0 0 -1 1 0 0 0 0 -1 1 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p + beta_dp - beta_f - beta_df > 0
    fprintf('F < P, p = %e\n', p);
else
    fprintf('F > P, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(6, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];






%% MCI

fprintf('\n----------------------------- MCI -----------------------------\n');

fprintf('\nBaseline\n\n');

% Overall: beta_s + beta_ms = beta_f + beta_mf = beta_p + beta_mp = 0
H = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0];
c = [0; 0; 0];
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall: p = %e\n', p);

% S - T: beta_s + beta_ms = 0
H = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_s + beta_ms > 0
    fprintf('T < S, p = %e\n', p);
else
    fprintf('T > S, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(7, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% F - T: beta_f + beta_mf = 0
H = [0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_f + beta_mf > 0
    fprintf('T < F, p = %e\n', p);
else
    fprintf('T > F, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(8, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% P - T: beta_p + beta_mp = 0
H = [0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p + beta_mp > 0
    fprintf('T < P, p = %e\n', p);
else
    fprintf('T > P, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(9, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% F - S:  beta_f + beta_mf - beta_s - beta_ms = 0
H = [0 0 0 -1 1 0 -1 1 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_f + beta_mf - beta_s - beta_ms > 0
    fprintf('S < F, p = %e\n', p);
else
    fprintf('S > F, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(10, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];

% P - S: beta_p + beta_mp - beta_s - beta_ms = 0
H = [0 0 0 -1 0 1 -1 0 1 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p + beta_mp - beta_s - beta_ms > 0
    fprintf('S < P, p = %e\n', p);
else
    fprintf('S > P, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(11, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];

% P - F: beta_p + beta_mp - beta_f - beta_mf = 0
H = [0 0 0 0 -1 1 0 -1 1 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p + beta_mp - beta_f - beta_mf > 0
    fprintf('F < P, p = %e\n', p);
else
    fprintf('F > P, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(12, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];






%% Normal

fprintf('\n----------------------------- Normal -----------------------------\n');

fprintf('\nBaseline\n\n');

% Overall: beta_s = beta_f = beta_p = 0
H = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
c = [0; 0; 0];
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall: p = %e\n', p);

% S - T: beta_s = 0
H = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_s > 0
    fprintf('T < S, p = %e\n', p);
else
    fprintf('T > S, p = %e\n', p);
end
idx1 = find(H==1);
baseline(13, :) = [betas(idx1), sqrt(stats.covb(idx1, idx1)), p];

% F - T: beta_f = 0
H = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_f > 0
    fprintf('T < F, p = %e\n', p);
else
    fprintf('T > F, p = %e\n', p);
end
idx1 = find(H==1);
baseline(14, :) = [betas(idx1), sqrt(stats.covb(idx1, idx1)), p];

% P - T: beta_p = 0
H = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p > 0
    fprintf('T < P, p = %e\n', p);
else
    fprintf('T > P, p = %e\n', p);
end
idx1 = find(H==1);
baseline(15, :) = [betas(idx1), sqrt(stats.covb(idx1, idx1)), p];

% F - S: beta_f - beta_s = 0
H = [0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_f - beta_s > 0
    fprintf('S < F, p = %e\n', p);
else
    fprintf('S > F, p = %e\n', p);
end
idx1 = find(H==1);
idx2 = find(H==-1);
baseline(16, :) = [betas(idx1)-betas(idx2), sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2)), p];

% P - S: beta_p - beta_s = 0
H = [0 0 0 -1 0 1 0 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p - beta_s > 0
    fprintf('S < P, p = %e\n', p);
else
    fprintf('S > P, p = %e\n', p);
end
idx1 = find(H==1);
idx2 = find(H==-1);
baseline(17, :) = [betas(idx1)-betas(idx2), sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2)), p];

% P - F: beta_p - beta_f = 0
H = [0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0];
c = 0;
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0 % beta_p - beta_f > 0
    fprintf('F < P, p = %e\n', p);
else
    fprintf('F > P, p = %e\n', p);
end
idx1 = find(H==1);
idx2 = find(H==-1);
baseline(18, :) = [betas(idx1)-betas(idx2), sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2)), p];
