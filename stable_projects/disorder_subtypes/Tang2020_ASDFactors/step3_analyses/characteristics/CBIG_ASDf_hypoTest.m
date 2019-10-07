function [baseline, pVals] = CBIG_ASDf_hypoTest(betas, stats, numReg, k)
% [baseline, pVals] = CBIG_ASDf_hypoTest(betas, stats, numReg, k)
% 
% This function performs hypothesis test for GLM for two-, three- and
% four-factor estimates.
% 
% Input:
%     - betas:
%           betas output from MATLAB glmfit.m function
%     - stats:
%           stats output from MATLAB glmfit.m function
%     - numReg:
%           Integer. Number of regressors
%     - k:
%           Integer. Number of factors
% Output:
%     - baseline:
%           Statistics needed to plot GLM results using herrorbar, which consists of
%           estimated difference between factors, standard error and p-value.
%     - pVals:
%           Vx1 vector, where V is the number of pairwise comparisons.
%           E.g., if k=3, V is 3; if k=4, V is 6. This is the pairwise
%           p-values of the hypothesis tests.
%
% Example:
%       [baseline, pVals] = CBIG_ASDf_hypoTest(betas, stats, 13, 3)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

pVals = []; % all pairwise comparisons

%%% K=2
if k == 2   
    % factor1 - factor2
    H = [0 1 zeros(1,numReg)];
    c = 0; % H*b = c
    
    p = linhyptest(betas, stats.covb, c, H, stats.dfe);
    
    % output to console
    if H*betas > 0
        fprintf('factor1 > factor2, p = %e\n', p);
    else
        fprintf('factor1 < factor2, p = %e\n', p);
    end
    
    idx1 = find(H==1);
    baseline = [H*betas sqrt(stats.covb(idx1, idx1)) p];
    pVals = [pVals; p];

%%% K=3
elseif k==3
    % Overall
    % Here, factor3 is implicitly modelled
    H = [0 1 0 zeros(1,numReg); 0 0 1 zeros(1,numReg)];
    c = [0; 0]; % H*b = c
    p = linhyptest(betas, stats.covb, c, H, stats.dfe);
    
    % output to console
    fprintf('Overall p = %e\n', p);
    
    % factor1 - factor3
    H = [0 1 0 zeros(1,numReg)];
    c = 0; % H*b = c
    p = linhyptest(betas, stats.covb, c, H, stats.dfe);
    if H*betas > 0
        fprintf('factor1 > factor3, p = %e\n', p);
    else
        fprintf('factor1 < factor3, p = %e\n', p);
    end
    idx1 = find(H==1);
    baseline(1, :) = [H*betas sqrt(stats.covb(idx1, idx1)) p];
    pVals = [pVals; p];
    
    % factor2 - factor3
    H = [0 0 1 zeros(1,numReg)];
    c = 0; % H*b = c
    p = linhyptest(betas, stats.covb, c, H, stats.dfe);
    if H*betas > 0
        fprintf('factor2 > factor3, p = %e\n', p);
    else
        fprintf('factor2 < factor3, p = %e\n', p);
    end
    idx1 = find(H==1);
    baseline(2, :) = [H*betas sqrt(stats.covb(idx1, idx1)) p];
    pVals = [pVals; p];
    
    % factor2 - factor1
    H = [0 -1 1 zeros(1,numReg)];
    c = 0; % H*b = c
    p = linhyptest(betas, stats.covb, c, H, stats.dfe);
    if H*betas > 0
        fprintf('factor2 > factor1, p = %e\n', p);
    else
        fprintf('factor2 < factor1, p = %e\n', p);
    end
    idx1 = find(H==1);
    idx2 = find(H==-1);
    baseline(3, :) = [H*betas...
        sqrt(stats.covb(idx1, idx1)+stats.covb(idx2, idx2)-2*stats.covb(idx1, idx2))...
        p];
    pVals = [pVals; p];
else
    fprintf('NOT CONFIGURED! Choose from K = 2 or K = 3.');
end

end
