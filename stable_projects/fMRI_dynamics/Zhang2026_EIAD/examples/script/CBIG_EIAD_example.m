function CBIG_EIAD_example(run_check)
% CBIG_EIAD_EXAMPLE  Quick check that the E/I-AD analysis code is set up correctly.
%
% Runs a reduced version of the ADNI analysis in CBIG_EIAD_ADNI.m on the
% shared participant table:
%   (1) GLM between the mean cortical E/I ratio and CSF Abeta42;
%   (2) region-wise GLM of the E/I ratio vs CSF Abeta42, and its alignment with
%       the sensorimotor-association (SA) axis;
%   (3) mediation of the CSF Abeta42 -> memory relationship by the mean cortical
%       E/I ratio (nonparametric bootstrap).
% The regional mediation and the figures of the full analysis are skipped, which
% is what keeps the run time short.
%
% NOTE ON THE SHARED DATA. To comply with the ADNI data-use agreement, the only
% real values released here are the subject IDs and the 68 harmonized regional
% E/I ratios. The AB42_CSF and PHCMemory columns in ADNI_df.csv are SYNTHETIC
% and are provided only so that this example runs end to end. The numbers this 
% script prints are therefore illustrative and are NOT the values reported in the paper. 
% To reproduce the published results, obtain the real AB42_CSF, PHCMemory, and covariates from
% ADNI using the released subject IDs (see ../README.md).
%
% Inputs (read from the current folder):
%   ADNI_df.csv          - participant table (subject ID, synthetic CSF Abeta42,
%                          synthetic PHC memory score, and 68 harmonized regional
%                          E/I ratios)
%   SA_rank_desikan.txt  - SA-axis rank of each of the 68 Desikan regions
%
% Input (optional):
%   run_check - if true (default), the results are compared against the reference
%               output at the end of the run. Pass false to only produce and save
%               the results without comparing (used by the unit test, which calls
%               CBIG_EIAD_check_example_results itself).
%
% Output:
%   ../output/example_output.mat - struct with the results of this run, compared
%                                  against ../reference_output/reference_output.mat
%
% When run_check is true, the results are compared field by field against that
% reference with a tolerance of 1e-6, by the helper function
% CBIG_EIAD_check_example_results.
%
% Dependencies (in this folder):
%   bootstrap_mediation             - local function (below)
%   CBIG_EIAD_check_example_results - compares the results to the reference output
%
% Requires the MATLAB Statistics and Machine Learning Toolbox (fitglm, zscore,
% randsample). Typical run time is under 10 seconds.
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 1 || isempty(run_check)
    run_check = true;
end

clc
close all

%% Load data
df         = readtable('ADNI_df.csv');
phc_memory = df.PHCMemory;   % synthetic
abeta42    = df.AB42_CSF;    % synthetic

% 68 harmonized regional E/I ratios (Desikan atlas)
EI = df{:, "HarmonizedEI" + string(1:68)};

% exclude participants with a missing CSF marker
idx = isnan(abeta42);
EI(idx, :)         = [];
abeta42(idx, :)    = [];
phc_memory(idx, :) = [];

output_struct = struct();

%% Analysis 1: mean-cortex GLM
EI_mean_z = zscore(mean(EI, 2));
abeta_z   = zscore(log10(abeta42));
mdl     = fitglm(abeta_z, EI_mean_z, 'Distribution', 'normal', 'Link', 'identity');
coefTbl = mdl.Coefficients;
rn = coefTbl.Properties.RowNames;
mapFrom = {'x1'};
mapTo   = {'log10(abeta42)'};
[lia, loc] = ismember(rn, mapFrom);
rn(lia) = mapTo(loc(lia));
coefTbl.Properties.RowNames = rn;
disp(coefTbl)
output_struct.GLM = coefTbl;

%% Analysis 2: regional GLM and SA-axis alignment
slope = zeros(68, 1);
for i = 1:68
    mdl = fitglm(abeta_z, zscore(EI(:, i)), 'Distribution', 'normal', 'Link', 'identity');
    slope(i) = mdl.Coefficients.Estimate(2);
end
output_struct.regional_slope = slope;

SA = dlmread('SA_rank_desikan.txt');
output_struct.regional_slope_SA_corr = corr(slope, SA, 'type', 'Spearman');

%% Analysis 3: mediation (CSF Abeta42 -> E/I ratio -> memory)
% restrict to participants with a PHC memory score
keep = ~isnan(phc_memory);
X = log10(abeta42(keep));       % independent variable
M = mean(EI(keep, :), 2);       % mediator
Y = phc_memory(keep);           % dependent variable

med = bootstrap_mediation(zscore(X), zscore(M), zscore(Y), []);
disp(['a = ' num2str(med.ahat) ', p = ' num2str(med.p_a)]);
disp(['b = ' num2str(med.bhat) ', p = ' num2str(med.p_b)]);
disp(['direct effect c prime = ' num2str(med.c_prime_hat) ', p = ' num2str(med.p_c_prime)]);
disp(['indirect effect ab = ' num2str(med.indirect_hat) ', p = ' num2str(med.p_ab)]);
disp(['proportion: ' num2str(med.indirect_hat/med.c_prime_hat)])
output_struct.global_a       = med.ahat;
output_struct.global_b       = med.bhat;
output_struct.global_c_prime = med.c_prime_hat;
output_struct.global_ab      = med.indirect_hat;

%% Save
if ~exist('../output', 'dir')
    mkdir('../output')
end
save('../output/example_output.mat', 'output_struct')

%% Compare against the reference output
if run_check
    CBIG_EIAD_check_example_results(output_struct)
end

end

%% ---------------------------------------------------------------------------
function result = bootstrap_mediation(X, M, Y, Z, B)
% BOOTSTRAP_MEDIATION  Nonparametric bootstrap test of mediation effects.
%
% X, M, Y are n x 1 vectors; Z is an n x p covariate matrix ([] if none);
% B is the number of bootstrap draws (default 10000).
%
% Returns a struct with the point estimates (ahat, bhat, c_hat, c_prime_hat,
% indirect_hat), their bootstrap distributions, two-sided bootstrap p-values,
% and the 95% percentile CI of the indirect effect.
%
% This is the same routine as in CBIG_EIAD_ADNI.m, repeated here so that the
% example runs standalone.

if nargin < 5, B = 10000; end
rng(1);  % reproducibility

n = numel(X);
Zmat = Z;   % covariate matrix
if isempty(Zmat)
    Zmat = zeros(n, 0);   % n-by-0 so row indexing (Zmat(idx,:)) still works
end

Xcol = X(:); Mcol = M(:); Ycol = Y(:);

% --- Fit on the original sample ---
% Path a: M ~ X + Z
coefA = regress(Mcol, [ones(n,1), Xcol, Zmat]);
ahat  = coefA(2);

% Path b and direct effect c': Y ~ M + X + Z
coefB       = regress(Ycol, [ones(n,1), Mcol, Xcol, Zmat]);
bhat        = coefB(2);
c_prime_hat = coefB(3);

% Total effect c: Y ~ X + Z
coefC = regress(Ycol, [ones(n,1), Xcol, Zmat]);
c_hat = coefC(2);

indirect_hat = ahat * bhat;

% --- Bootstrap ---
boot_a      = nan(B,1);
boot_b      = nan(B,1);
boot_cprime = nan(B,1);
boot_c      = nan(B,1);
boot_ab     = nan(B,1);

for b = 1:B
    idx = randsample(n, n, true);
    coefA_b = regress(Mcol(idx), [ones(n,1), Xcol(idx), Zmat(idx,:)]);
    coefB_b = regress(Ycol(idx), [ones(n,1), Mcol(idx), Xcol(idx), Zmat(idx,:)]);
    coefC_b = regress(Ycol(idx), [ones(n,1), Xcol(idx), Zmat(idx,:)]);

    boot_a(b)      = coefA_b(2);
    boot_b(b)      = coefB_b(2);
    boot_cprime(b) = coefB_b(3);
    boot_c(b)      = coefC_b(2);
    boot_ab(b)     = coefA_b(2) * coefB_b(2);
end

% --- Two-sided bootstrap p-values ---
pval_fun = @(bootDist, est) 2 * min(mean(bootDist >= est), mean(bootDist <= est));

result.ahat         = ahat;
result.bhat         = bhat;
result.c_hat        = c_hat;
result.c_prime_hat  = c_prime_hat;
result.indirect_hat = indirect_hat;

result.p_a       = pval_fun(boot_a, 0);
result.p_b       = pval_fun(boot_b, 0);
result.p_c       = pval_fun(boot_c, 0);
result.p_c_prime = pval_fun(boot_cprime, 0);
result.p_ab      = pval_fun(boot_ab, 0);

result.boot_a      = boot_a;
result.boot_b      = boot_b;
result.boot_cprime = boot_cprime;
result.boot_c      = boot_c;
result.boot_ab     = boot_ab;

% 95% percentile CI for the indirect effect
alpha = 0.05;
result.CI_ab = prctile(boot_ab, [100*alpha/2, 100*(1-alpha/2)]);
result.B = B;

end
