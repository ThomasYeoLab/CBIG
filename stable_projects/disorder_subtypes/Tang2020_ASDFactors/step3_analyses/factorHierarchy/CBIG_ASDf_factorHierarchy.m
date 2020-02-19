function [corr_brainMap, corr_factorComp] = CBIG_ASDf_factorHierarchy(k, inputDir, r_k, r_kplus)
% [corr_brainMap, corr_factorComp] = CBIG_ASDf_factorHierarchy(k, inputDir, r_k, r_kplus)
% 
% This function quantifies the nested hierarchy of factors in terms of both 
% E(RSFC patterns|Factor) and Pr(Factor|Participant). 
% E(RSFC patterns|Factor) is the factor-specific hypo/hyper RSFC patterns, 
% and is computed by beta*(2*rho-1).
% 
% This function will explore all possible combinations of splits, and
% output the quantifications on the console.
% 
% NOTE: the factors here are not re-ordered. Therefore, in the text displayed
% on console, 'Factor1' means the 1-st factor in the specified run (random initialization),
% and does not necessarily match factor 1 in our paper.
%
% Input:
%     - k:
%           Integer. Number of factors from which one factor splits
%     - inputDir:
%           Full path to the directory where the model estimation results
%           are located
%     - r_k:
%           Integer. The run number for k-factor model
%     - r_kplus:
%           Integer. The run number for (k+1)-factor model
%
% Output:
%     - corr_brainMap:
%           Nx1 vector, where N is the number of possible combinations of splits. 
%           The averaged correlation for E(RSFC patterns|Factor) in all
%           possible combinations of splits.
%     - corr_factorComp:
%           Nx1 vector, where N is the number of possible combinations of splits.
%           The averaged correlation for Pr(Factor|Participant) in all possible 
%           combinations of splits. 
%
% Example:
%     [corr_brainMap, corr_factorComp] = CBIG_ASDf_factorHierarchy(2, '~/example_output/estimate', 96, 94)
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step2_polarLDA'));

%% Number of possible combinations of splits
num_pairs = nchoosek(k+1,2);
combinations = combnk(1:(k+1),2);

corr_brainMap = zeros(num_pairs,1);
corr_factorComp = zeros(num_pairs,1);

%% Compute E(FC patterns | Factor) for k-factor and (k+1)-factor models
beta_k = exp(load(fullfile(inputDir,sprintf('k%s/r%s', num2str(k), num2str(r_k)),'final.beta')));
rho_k = exp(load(fullfile(inputDir,sprintf('k%s/r%s', num2str(k), num2str(r_k)),'final.rho')));
Mean_k = beta_k.*(2*rho_k-1);

beta_kplus = exp(load(fullfile(inputDir,sprintf('k%s/r%s', num2str(k+1), num2str(r_kplus)),'final.beta')));
rho_kplus = exp(load(fullfile(inputDir,sprintf('k%s/r%s', num2str(k+1), num2str(r_kplus)),'final.rho')));
Mean_kplus = beta_kplus.*(2*rho_kplus-1);

%% Compute Pr(Factor | Participant) for k-factor and (k+1)-factor models
gamma_k = load(fullfile(inputDir,['k' num2str(k)],['r' num2str(r_k)],'final.gamma'));
factorComp_k = bsxfun(@times, gamma_k, 1./(sum(gamma_k, 2)));

gamma_kplus = load(fullfile(inputDir,['k' num2str(k+1)],['r' num2str(r_kplus)],'final.gamma'));
factorComp_kplus = bsxfun(@times, gamma_kplus, 1./(sum(gamma_kplus, 2)));

%% Now perform exhaustive search over all possible splits
idx_pair = 1;

for p = 1:size(combinations,1)
    Mean_comb = (Mean_kplus(combinations(p,1),:) + Mean_kplus(combinations(p,2),:))/2;
    factorComp_comb = (factorComp_kplus(:,combinations(p,1)) + factorComp_kplus(:,combinations(p,2)))/2;
    Mean_kplus_comb = zeros(size(Mean_k));
    factorComp_kplus_comb = zeros(size(factorComp_k));
    
    index = ~ismember(1:(k+1), combinations(p,:));
    
    Mean_kplus_comb(1:(k-1),:) = Mean_kplus(index,:); % keep other factors
    Mean_kplus_comb(k,:) = Mean_comb; % the combined factor is the last index in Mean_kplus_comb
    factorComp_kplus_comb(:,1:(k-1)) = factorComp_kplus(:,index);
    factorComp_kplus_comb(:,k) = factorComp_comb;
    
    order = CBIG_ASDf_hunMatch(k, Mean_k, Mean_kplus_comb);
    idx_splitted = find(order == k); % which factor splits?
    
    sumCorr_brainMap = 0;
    sumCorr_factorComp = 0;
    for idx = 1:k
        corrMat_brainMap = corrcoef(Mean_k(idx,:), Mean_kplus_comb(order(idx),:));
        corrMat_factorComp = corrcoef(factorComp_k(:,idx), factorComp_kplus_comb(:,order(idx)));
        
        sumCorr_brainMap = sumCorr_brainMap + corrMat_brainMap(1,2);
        sumCorr_factorComp = sumCorr_factorComp + corrMat_factorComp(1,2);
    end
    avgCorr_brainMap = sumCorr_brainMap/k;
    avgCorr_factorComp = sumCorr_factorComp/k;
    
    fprintf('Assuming Factor%d and Factor%d in %d-factor model are splits from Factor%d in %d-factor model...\n', ...
combinations(p,1),combinations(p,2),k+1,idx_splitted,k);
    fprintf('Averaged correlation of E(FC patterns|Factor) is: %f\n', avgCorr_brainMap);
    fprintf('Averaged correlation of Pr(Factor|Participant) is : %f\n\n', avgCorr_factorComp);
    
    corr_brainMap(idx_pair,:) = avgCorr_brainMap;
    corr_factorComp(idx_pair,:) = avgCorr_factorComp;
    
    idx_pair = idx_pair + 1;
end

%% Remove path
rmpath(fullfile(CODE_DIR,'step2_polarLDA'));
