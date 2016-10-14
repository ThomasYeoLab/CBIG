function baseline = CBIG_hypoTest_acrossStages(betas, stats)

% baseline = CBIG_hypoTest_acrossStages(betas, stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fprintf('\n\n****************** Combine Factors, Across Stages ******************\n\n');


baseline = zeros(3, 3);

% Overall: beta_m = beta_d = 0
H = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% MCI - CN: beta_m = 0
H = [0 1 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('MCI > CN, p = %e\n', p);
else
    fprintf('MCI < CN, p = %e\n', p);
end
idx = find(H==1);
baseline(1, :) = [betas(idx), sqrt(stats.covb(idx, idx)), p];

% AD - MCI: beta_d - beta_m = 0
H = [0 -1 1 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > MCI, p = %e\n', p);
else
    fprintf('AD < MCI, p = %e\n', p);
end
idx1 = find(H==1);
idx2 = find(H==-1);
baseline(2, :) = [...
    betas(idx1)-betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2)),...
    p];

% AD - CN: beta_d = 0
H = [0 0 1 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > CN, p = %e\n', p);
else
    fprintf('AD < CN, p = %e\n', p);
end
idx = find(H==1);
baseline(3, :) = [betas(idx), sqrt(stats.covb(idx, idx)), p];