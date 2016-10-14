function baseline = CBIG_hypoTest_3t_withinFactorAcrossStages(betas, stats)

% baseline = CBIG_hypoTest_3t_withinFactorAcrossStages(betas, stats)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fprintf('\n\n****************** Within Factor, Across Stages ******************\n\n');


baseline = zeros(9, 3);

%% Temporal

fprintf('\n--------------------------------- Temporal ---------------------------------\n\n');

% Overall: beta_m = beta_d = 0
H = [0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% MCI - CN: beta_m = 0
H = [0 1 0 0 0 0 0 0 0 0 0 0 0 0];
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
H = [0 -1 1 0 0 0 0 0 0 0 0 0 0 0];
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
H = [0 0 1 0 0 0 0 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > CN, p = %e\n', p);
else
    fprintf('AD < CN, p = %e\n', p);
end
idx = find(H==1);
baseline(3, :) = [betas(idx), sqrt(stats.covb(idx, idx)), p];


%% Subcortical

fprintf('\n--------------------------------- Subcortical ---------------------------------\n\n');

% Overall: beta_m + beta_ms = beta_d + beta_ds = 0
H = [0 1 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 1 0 0 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% MCI - CN: beta_m + beta_ms = 0
H = [0 1 0 0 0 1 0 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('MCI > CN, p = %e\n', p);
else
    fprintf('MCI < CN, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(4, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% AD - MCI: beta_d + beta_ds - beta_m - beta_ms = 0
H = [0 -1 1 0 0 -1 0 1 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > MCI, p = %e\n', p);
else
    fprintf('AD < MCI, p = %e\n', p);
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

% AD - CN: beta_d + beta_ds = 0
H = [0 0 1 0 0 0 0 1 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > CN, p = %e\n', p);
else
    fprintf('AD < CN, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(6, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

%% Cortical

fprintf('\n--------------------------------- Cortical ---------------------------------\n\n');

% Overall: beta_m + beta_mc = beta_d + beta_dc = 0
H = [0 1 0 0 0 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 1 0 0 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% MCI - CN: beta_m + beta_mc = 0
H = [0 1 0 0 0 0 1 0 0 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('MCI > CN, p = %e\n', p);
else
    fprintf('MCI < CN, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(7, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];

% AD - MCI: beta_d + beta_dc - beta_m - beta_mc = 0
H = [0 -1 1 0 0 0 -1 0 1 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > MCI, p = %e\n', p);
else
    fprintf('AD < MCI, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
ind = find(H==-1);
idx3 = ind(1);
idx4 = ind(2);
baseline(8, :) = [...
    betas(idx1)+betas(idx2)-betas(idx3)-betas(idx4),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+stats.covb(idx3, idx3)+stats.covb(idx4, idx4)+2*stats.covb(idx1, idx2)-2*stats.covb(idx1, idx3)-2*stats.covb(idx1, idx4)-2*stats.covb(idx2, idx3)-2*stats.covb(idx2, idx4)+2*stats.covb(idx3, idx4)),...
    p];

% AD - CN: beta_d + beta_dc = 0
H = [0 0 1 0 0 0 0 0 1 0 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('AD > CN, p = %e\n', p);
else
    fprintf('AD < CN, p = %e\n', p);
end
ind = find(H==1);
idx1 = ind(1);
idx2 = ind(2);
baseline(9, :) = [...
    betas(idx1)+betas(idx2),...
    sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)+2*stats.covb(idx1, idx2)),...
    p];