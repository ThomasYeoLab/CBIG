function baseline = CBIG_hypoTest_3t(betas, stats)

% baseline = CBIG_hypoTest_3t(betas, stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

baseline = zeros(9, 3);


%% AD

fprintf('\n--------------------------------- AD ---------------------------------\n\n');

% Overall: beta_s + beta_ds = beta_c + beta_dc = 0
H = [0 0 0 1 0 0 0 1 0 0 0 0 0; 0 0 0 0 1 0 0 0 1 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% S - T: beta_s + beta_ds = 0
H = [0 0 0 1 0 0 0 1 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('S > T, p = %e\n', p);
else
    fprintf('S < T, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(1, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% C - T: beta_c + beta_dc = 0
H = [0 0 0 0 1 0 0 0 1 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T, p = %e\n', p);
else
    fprintf('C < T, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(2, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% C - S: beta_c - beta_s + beta_dc - beta_ds = 0
H = [0 0 0 -1 1 0 0 -1 1 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > S, p = %e\n', p);
else
    fprintf('C < S, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(3, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];

%% MCI

fprintf('\n--------------------------------- MCI ---------------------------------\n\n');

% Overall: beta_s + beta_ms = beta_c + beta_mc = 0
H = [0 0 0 1 0 1 0 0 0 0 0 0 0; 0 0 0 0 1 0 1 0 0 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% S - T: beta_s + beta_ms = 0
H = [0 0 0 1 0 1 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('S > T, p = %e\n', p);
else
    fprintf('S < T, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(4, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% C - T: beta_c + beta_mc = 0
H = [0 0 0 0 1 0 1 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T, p = %e\n', p);
else
    fprintf('C < T, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(5, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% C - S: beta_c - beta_s + beta_mc - beta_ms = 0
H = [0 0 0 -1 1 -1 1 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > S, p = %e\n', p);
else
    fprintf('C < S, p = %e\n', p);
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

%% CN

fprintf('\n--------------------------------- CN ---------------------------------\n\n');

% Overall: beta_s = beta_c = 0
H = [0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% S - T: beta_s = 0
H = [0 0 0 1 0 0 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('S > T, p = %e\n', p);
else
    fprintf('S < T, p = %e\n', p);
end
idx1 = find(H==1);
baseline(7, :) = [betas(idx1), sqrt(stats.covb(idx1, idx1)), p];

% C - T: beta_c = 0
H = [0 0 0 0 1 0 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T, p = %e\n', p);
else
    fprintf('C < T, p = %e\n', p);
end
idx1 = find(H==1);
baseline(8, :) = [betas(idx1), sqrt(stats.covb(idx1, idx1)), p];

% C - S: beta_c - beta_s = 0
H = [0 0 0 -1 1 0 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > S, p = %e\n', p);
else
    fprintf('C < S, p = %e\n', p);
end
idx1 = find(H==1);
idx2 = find(H==-1);
baseline(9, :) = [...
    betas(idx1)-betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2)),...
    p];
