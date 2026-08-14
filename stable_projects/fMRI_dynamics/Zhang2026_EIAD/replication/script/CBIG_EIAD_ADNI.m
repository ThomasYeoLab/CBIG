function CBIG_EIAD_ADNI()
% CBIG_EIAD_ADNI  Reproduce the main ADNI analyses of the E/I-AD study.
%
% On the ADNI cohort, this wrapper runs:
%   (1) GLM between the mean cortical E/I ratio and CSF Abeta42;
%   (2) region-wise GLM of the E/I ratio vs CSF Abeta42, and its alignment with
%       the sensorimotor-association (SA) axis;
%   (3) mediation of the CSF Abeta42 -> cognition relationship by the E/I ratio,
%       at the mean-cortex and regional levels (nonparametric bootstrap).
%
% Inputs (read from the current folder):
%   ADNI_df.csv          - participant table (demographics, CSF Abeta42, PHC
%                          memory score, and 68 harmonized regional E/I ratios)
%   SA_rank_desikan.txt  - SA-axis rank of each of the 68 Desikan regions
%
% Output:
%   ../output/ADNI.mat   - struct with the results reported in the paper.
%                          A reference copy is under ../reference_output/.
%
% Dependencies:
%   bootstrap_mediation        - local function (below)
%   CBIG_EIAD_draw_surface_map - only for the (commented) surface figures
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

clc
close all

%% Load data
df         = readtable('ADNI_df.csv');

% the released table ships with the covariate columns empty; stop with a clear
% message if any of them has not been repopulated (see ../README.md)
covVars = {'age', 'sex_F0M1', 'FD', 'yrsOfEducation'};
emptyCol = varfun(@(v) all(ismissing(v)), df(:, covVars), 'OutputFormat', 'uniform');
if any(emptyCol)
    error(['One or more covariate columns (%s) are empty. The released ' ...
           'ADNI_df.csv ships without demographic, biomarker, or cognitive ' ...
           'values; repopulate every covariate column from ADNI using the ' ...
           'subject IDs before running (see the "Repopulating the tables" ' ...
           'section of ../README.md).'], strjoin(covVars(emptyCol), ', '))
end

age        = df.age;
sex        = df.sex_F0M1;
FD         = df.FD;
education  = df.yrsOfEducation;
phc_memory = df.PHCMemory;
abeta42    = df.AB42_CSF;

% 68 harmonized regional E/I ratios (Desikan atlas)
EI  = df{:, "HarmonizedEI" + string(1:68)};
cov = [age sex FD education];

% exclude participants with any missing covariate or CSF marker
idx = any(isnan([cov abeta42]), 2);
EI(idx, :)      = [];
cov(idx, :)     = [];
abeta42(idx, :) = [];

% continuous covariates (age, FD, education) are z-scored; sex (2nd entry) is not
cont_flag = [1 0 1 1];

output_struct = struct();

%% Analysis 1: mean-cortex GLM
EI_mean = mean(EI, 2);
abeta   = log10(abeta42);

% E/I vs CSF Abeta42, adjusting for age, sex, FD and education
EI_residual    = residualize(EI_mean, cov);
abeta_residual = residualize(abeta,   cov);
[r, p] = corr(EI_residual, abeta_residual);
disp(['E/I - CSF correlation: r = ' num2str(r) ', p = ' num2str(p)])
output_struct.EI_CSF_corr = r;
figure(1)
local_scatter_fit(abeta_residual, EI_residual, ...
    'log10(CSF Abeta42) residual', 'E/I ratio residual')

% E/I vs age, adjusting for CSF Abeta42, sex, FD and education
EI_residual  = residualize(EI_mean,   [abeta cov(:, 2:end)]);
age_residual = residualize(cov(:, 1), [abeta cov(:, 2:end)]);
[r, p] = corr(EI_residual, age_residual);
disp(['E/I - age correlation: r = ' num2str(r) ', p = ' num2str(p)])
output_struct.EI_age_corr = r;
figure(2)
local_scatter_fit(age_residual, EI_residual, 'age residual', 'E/I ratio residual')

% GLM on the (non-residualized) mean E/I ratio while controlling for covariates
cov_z     = zscore_continuous(cov, cont_flag);
EI_mean_z = zscore(mean(EI, 2));
abeta_z   = zscore(log10(abeta42));
mdl       = fitglm([abeta_z cov_z], EI_mean_z, 'Distribution', 'normal', 'Link', 'identity');
coefTbl   = mdl.Coefficients;
rn = coefTbl.Properties.RowNames;
mapFrom = {'x1','x2','x3','x4','x5'};
mapTo   = {'log10(abeta42)','age','sex','FD','yrs education'};
[lia, loc] = ismember(rn, mapFrom);
rn(lia) = mapTo(loc(lia));
coefTbl.Properties.RowNames = rn;
disp(coefTbl)
output_struct.GLM = coefTbl;

%% Analysis 2: regional GLM and SA-axis alignment
slope   = zeros(68, 1);
slope_p = zeros(68, 1);
for i = 1:68
    mdl = fitglm([abeta_z cov_z], zscore(EI(:, i)), 'Distribution', 'normal', 'Link', 'identity');
    slope(i)   = mdl.Coefficients.Estimate(2);
    slope_p(i) = mdl.Coefficients.pValue(2);
end
output_struct.regional_slope = slope;

% surface visualization (requires CBIG_EIAD_draw_surface_map)
roi_list = ones(72, 1);
roi_list([1 5 37 41]) = 0;
% CBIG_EIAD_draw_surface_map(slope, roi_list, 'regional slope')

% alignment of regional slopes with the SA axis
SA = dlmread('SA_rank_desikan.txt');
figure(4)
local_scatter_fit(SA, slope, 'SA ranking', 'slope')
output_struct.regional_slope_SA_corr = corr(slope, SA, 'type', 'Spearman');

%% Analysis 3: mediation (CSF Abeta42 -> E/I ratio -> cognition)
% restrict to participants with a PHC memory score
keep = ~isnan(phc_memory);
X   = log10(abeta42(keep));                        % independent variable
Y   = phc_memory(keep);                            % dependent variable
Z_z = zscore_continuous(cov(keep, :), cont_flag);  % covariates

% mean-cortex mediation
M   = mean(EI(keep, :), 2);                        % mediator
med = bootstrap_mediation(zscore(X), zscore(M), zscore(Y), Z_z);
disp(['a = ' num2str(med.ahat) ', p = ' num2str(med.p_a)]);
disp(['b = ' num2str(med.bhat) ', p = ' num2str(med.p_b)]);
disp(['direct effect c prime = ' num2str(med.c_prime_hat) ', p = ' num2str(med.p_c_prime)]);
disp(['indirect effect ab = ' num2str(med.indirect_hat) ', p = ' num2str(med.p_ab)]);
disp(['proportion: ' num2str(med.indirect_hat/med.c_prime_hat)])
output_struct.global_a       = med.ahat;
output_struct.global_b       = med.bhat;
output_struct.global_c_prime = med.c_prime_hat;
output_struct.global_ab      = med.indirect_hat;
close all

% regional mediation
b    = zeros(68, 1);
b_p  = zeros(68, 1);
ab   = zeros(68, 1);
ab_p = zeros(68, 1);
for i = 1:68
    M   = EI(keep, i);
    med = bootstrap_mediation(zscore(X), zscore(M), zscore(Y), Z_z);
    b(i)    = med.bhat;
    b_p(i)  = med.p_b;
    ab(i)   = med.indirect_hat;
    ab_p(i) = med.p_ab;
end
output_struct.regional_b  = b;
output_struct.regional_ab = ab;
close all

%% Save
if ~exist('../output', 'dir')
    mkdir('../output')
end
save('../output/ADNI.mat', 'output_struct')

% surface visualization of regional effects, thresholded by FDR
% (requires CBIG_EIAD_draw_surface_map)
% p_adj = mafdr(b_p, 'BHFDR', 1);
% roi_list = p_adj < 0.05;
% b = b(roi_list);
% roi_list = [0; roi_list(1:3); 0; roi_list(4:34); 0; roi_list(35:37); 0; roi_list(38:68)];
% CBIG_EIAD_draw_surface_map(b, roi_list, 'b')
% p_adj = mafdr(ab_p, 'BHFDR', 1);
% roi_list = p_adj < 0.05;
% ab = ab(roi_list);
% roi_list = [0; roi_list(1:3); 0; roi_list(4:34); 0; roi_list(35:37); 0; roi_list(38:68)];
% CBIG_EIAD_draw_surface_map(ab, roi_list, 'ab')

end

%% ---------------------------------------------------------------------------
function res = residualize(y, X)
% Residuals of y after regressing out an intercept and the columns of X.
X   = [ones(size(X, 1), 1), X];
res = y - X * regress(y, X);
end

%% ---------------------------------------------------------------------------
function Xz = zscore_continuous(X, cont_flag)
% Z-score the continuous columns of X (cont_flag == 1); leave the rest unchanged.
Xz = X;
for i = 1:numel(cont_flag)
    if cont_flag(i) == 1
        Xz(:, i) = zscore(X(:, i));
    end
end
end

%% ---------------------------------------------------------------------------
function local_scatter_fit(x, y, xlab, ylab)
% Scatter of y vs x with a least-squares fit line (as used for the paper figures).
xlabel(xlab)
ylabel(ylab)
hold on
pfit  = polyfit(x, y, 1);
x_fit = linspace(min(x), max(x), 100);
y_fit = polyval(pfit, x_fit);
plot(x, y, '.', 'Color', [105 105 105]/255, 'MarkerSize', 12)
set(gca, 'LineWidth', 2)
set(gca, 'box', 'off')
plot(x_fit, y_fit, 'r-', 'LineWidth', 3)
hold off
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

if nargin < 5, B = 10000; end
rng(1);  % reproducibility

n = numel(X);
Zmat = Z;   % [] if no covariates

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
