function CBIG_hypoTest_3t(betas, stats)

% CBIG_hypoTest_3t(betas, stats)
% Overall
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

H = [0 1 0 0 0 0; 0 0 1 0 0 0];
c = [0; 0]; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
fprintf('Overall p = %e\n', p);

% S - T 
H = [0 1 0 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('S > T, p = %e\n', p);
else
    fprintf('S < T, p = %e\n', p);
end

% C - T
H = [0 0 1 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > T, p = %e\n', p);
else
    fprintf('C < T, p = %e\n', p);
end

% C - S
H = [0 -1 1 0 0 0];
c = 0; % H*b = c
p = linhyptest(betas, stats.covb, c, H, stats.dfe);
if H*betas > 0
    fprintf('C > S, p = %e\n', p);
else
    fprintf('C < S, p = %e\n', p);
end