function CBIG_fitGLM_hypoTest(rid_prob, rid_icv_vol, subCol, corCol)

% CBIG_fitGLM_hypoTest(rid_prob, rid_icv_vol, subCol, corCol)
%--------------------------------------------- Fit GLM
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% Construct X and y
noSubs = size(rid_icv_vol, 1);
X = zeros(noSubs, 3);
for idx = 1:noSubs
    % Logical indices for 810 subjects
    rowIdx = rid_prob(:, 1)==rid_icv_vol(idx, 1);
    s = rid_prob(rowIdx, subCol);
    c = rid_prob(rowIdx, corCol);
    i = rid_icv_vol(idx, 2);
    % X
    X(idx, :) = [s c i];
end
% y
y = rid_icv_vol(:, 3);

% Fit
[betas, ~, stats] = glmfit(X, y);

%--------------------------------------------- Hypothesis test

fprintf('\n');

% S - T
sCoef = betas(2);
if sCoef < 0
    fprintf('S < T\n');
else
    fprintf('S > T\n');
end
fprintf('p = %d\n', stats.p(2));

% C - T
cCoef = betas(3);
if cCoef < 0
    fprintf('C < T\n');
else
    fprintf('C > T\n');
end
fprintf('p = %d\n', stats.p(3));

% C - S
H = [0 -1 1 0];
c = 0;
if H*betas > 0
    fprintf('C > S\n');
else
    fprintf('C < S\n');
end
fprintf('p = %d\n', linhyptest(betas, stats.covb, c, H, stats.dfe));
