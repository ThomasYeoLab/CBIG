clear;
clc;
close all;

addpath(genpath('../functions/'));

CBIG_plotSetup;

K_0 = 2;
K = 10;

k_bestRun = cell(K-K_0+1, 2);
idx = 1;
for k = K_0:K
    r = CBIG_find_maxLikeRun(['../../lda/outputs/188blAD' sprintf('/k%i/', k)]);
    k_bestRun{idx, 1} = k;
    k_bestRun{idx, 2} = ['../../lda/outputs/188blAD/' sprintf('k%i/r%s/', k, r)];
    idx = idx+1;
end

CBIG_factorHierachy(k_bestRun);
