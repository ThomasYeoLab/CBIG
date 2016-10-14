function baseline = CBIG_hypoTest_2t(betas, stats)

% baseline = CBIG_hypoTest_2t(betas, stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

baseline = zeros(3, 3);


%% AD

fprintf('\n--------------------------------- AD ---------------------------------\n\n');

% C - T+S: beta_c + beta_dc = 0
H = [0 0 0 1 0 1 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T+S, p = %e\n', p);
else
    fprintf('C < T+S, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(1, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

%% MCI

fprintf('\n--------------------------------- MCI ---------------------------------\n\n');

% C - T+S: beta_c + beta_mc = 0
H = [0 0 0 1 1 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T+S, p = %e\n', p);
else
    fprintf('C < T+S, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(2, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

%% CN

fprintf('\n--------------------------------- CN ---------------------------------\n\n');

% C - T+S: beta_c = 0
H = [0 0 0 1 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T+S, p = %e\n', p);
else
    fprintf('C < T+S, p = %e\n', p);
end
idx1 = find(H==1);
baseline(3, :) = [betas(idx1), sqrt(stats.covb(idx1, idx1)), p];
