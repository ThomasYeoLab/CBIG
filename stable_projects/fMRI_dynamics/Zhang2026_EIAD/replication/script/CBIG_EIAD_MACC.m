function CBIG_EIAD_MACC()
% CBIG_EIAD_MACC  Reproduce the main MACC analyses of the E/I-AD study.
%
% On the MACC cohort, this wrapper runs:
%   (1) GLM between the mean cortical E/I ratio and plasma p-tau217;
%   (2) region-wise GLM of the E/I ratio vs plasma p-tau217, and its alignment
%       with the sensorimotor-association (SA) axis;
%   (3) mediation of the plasma p-tau217 -> cognition relationship by the E/I
%       ratio, at the mean-cortex and regional levels (nonparametric bootstrap).
%       Cognition is a single-factor score over the verbal-memory and MMSE items.
%
% Inputs (read from the current folder):
%   MACC_df.csv          - participant table (demographics, plasma p-tau217,
%                          verbal-memory and MMSE items, the 68 regional E/I
%                          ratios, and the underlying S_E and S_I means)
%   SA_rank_desikan.txt  - SA-axis rank of each of the 68 Desikan regions
%
% Output:
%   ../output/MACC.mat   - struct with the results reported in the paper.
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
df        = readtable('MACC_df.csv');

% the released table ships with the covariate columns empty; stop with a clear
% message if any of them has not been repopulated (see ../README.md)
covVars = {'age', 'sex_F0M1', 'FD', 'yrsOfEducation'};
emptyCol = varfun(@(v) all(ismissing(v)), df(:, covVars), 'OutputFormat', 'uniform');
if any(emptyCol)
    error(['One or more covariate columns (%s) are empty. The released ' ...
           'MACC_df.csv ships without demographic, biomarker, or cognitive ' ...
           'values; repopulate every covariate column from MACC using the ' ...
           'subject IDs before running (see the "Repopulating the tables" ' ...
           'section of ../README.md).'], strjoin(covVars(emptyCol), ', '))
end

age       = df.age;
sex       = df.sex_F0M1;
FD        = df.FD;
education = df.yrsOfEducation;
ptau217   = df.ptau217;

% 68 regional E/I ratios (Desikan atlas)
EI        = df{:, "EI" + string(1:68)};
biomarker = log10(ptau217);
cov       = [age sex FD education];

% exclude participants with any missing covariate or plasma marker
idx = any(isnan([cov biomarker]), 2);
EI(idx, :)        = [];
cov(idx, :)       = [];
biomarker(idx, :) = [];

% continuous covariates (age, FD, education) are z-scored; sex (2nd entry) is not
cont_flag = [1 0 1 1];

output_struct = struct();

%% Analysis 1: mean-cortex GLM
EI_mean = mean(EI, 2);

% E/I vs plasma p-tau217, adjusting for age, sex, FD and education
EI_residual        = residualize(EI_mean,   cov);
biomarker_residual = residualize(biomarker, cov);
[r, p] = corr(EI_residual, biomarker_residual);
disp(['E/I - plasma correlation: r = ' num2str(r) ', p = ' num2str(p)])
output_struct.EI_plasma_corr = r;
figure(1)
local_scatter_fit(biomarker_residual, EI_residual, ...
    'log10(ptau217) residual', 'E/I ratio residual')

% E/I vs head motion, adjusting for p-tau217, age, sex and education
EI_residual = residualize(EI_mean,   [biomarker cov(:, [1 2 4])]);
FD_residual = residualize(cov(:, 3), [biomarker cov(:, [1 2 4])]);
[r, p] = corr(EI_residual, FD_residual);
disp(['E/I - FD correlation: r = ' num2str(r) ', p = ' num2str(p)])
output_struct.EI_FD_corr = r;
figure(2)
local_scatter_fit(FD_residual, EI_residual, 'FD residual', 'E/I ratio residual')

% GLM on the (non-residualized) mean E/I ratio while controlling for covariates
cov_z       = zscore_continuous(cov, cont_flag);
EI_mean_z   = zscore(mean(EI, 2));
biomarker_z = zscore(biomarker);
mdl     = fitglm([biomarker_z cov_z], EI_mean_z, 'Distribution', 'normal', 'Link', 'identity');
coefTbl = mdl.Coefficients;
rn = coefTbl.Properties.RowNames;
mapFrom = {'x1','x2','x3','x4','x5'};
mapTo   = {'log10(ptau217)','age','sex','FD','yrs education'};
[lia, loc] = ismember(rn, mapFrom);
rn(lia) = mapTo(loc(lia));
coefTbl.Properties.RowNames = rn;
disp(coefTbl)
output_struct.GLM = coefTbl;

%% Analysis 2: regional GLM and SA-axis alignment
s   = zeros(68, 1);
s_p = zeros(68, 1);
for i = 1:68
    mdl = fitglm([biomarker_z cov_z], zscore(EI(:, i)), 'Distribution', 'normal', 'Link', 'identity');
    s(i)   = mdl.Coefficients.Estimate(2);
    s_p(i) = mdl.Coefficients.pValue(2);
end
output_struct.regional_slope = s;

% surface visualization (requires CBIG_EIAD_draw_surface_map)
roi_list = ones(72, 1);
roi_list([1 5 37 41]) = 0;
% CBIG_EIAD_draw_surface_map(s, roi_list, 'regional slope')

% alignment of regional slopes with the SA axis
SA = dlmread('SA_rank_desikan.txt');
figure(4)
local_scatter_fit(SA, s, 'SA ranking', 'slope')
output_struct.regional_slope_SA_corr = corr(s, SA, 'type', 'Spearman');

%% Analysis 3: mediation (plasma p-tau217 -> E/I ratio -> cognition)
% memory items of the participants retained above (same rows kept in idx)
mask = ~idx;

% word-list and story recall: immediate and delayed are summed within task
verbal_memory = df{mask, {'vbmwlrimdtrecal', 'vbmwlrdelyrecal', ...
                          'vbmsraimdt',      'vbmsradelyrecal'}};
verbal_memory = [sum(verbal_memory(:, [1 2]), 2), sum(verbal_memory(:, [3 4]), 2)];

% MMSE: the 10 orientation items are summed, then registration + recall
MMSE = df{mask, {'orientday_BL', 'orienttdydate_BL', 'orientcurmonth_BL', ...
                 'orientcuryear_BL', 'orientplaceclinic_BL', 'orienttime_BL', ...
                 'orientnameplace_BL', 'orientestate_BL', 'orientfloor_BL', ...
                 'orientcountry_BL', 'registrationscore_BL', 'recall_BL'}};
MMSE = [sum(MMSE(:, 1:10), 2), sum(MMSE(:, 11:12), 2)];

% single-factor score across the four memory subscores
memory_score = [verbal_memory MMSE];
[~, ~, ~, ~, F] = factoran(zscore(memory_score), 1, 'rotate', 'promax', 'scores', 'regression');
latent_memory = F;

X   = biomarker;                          % independent variable
Y   = latent_memory;                      % dependent variable
Z_z = zscore_continuous(cov, cont_flag);  % covariates

% mean-cortex mediation
M   = mean(EI, 2);                        % mediator
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

% regional mediation
b    = zeros(68, 1);
b_p  = zeros(68, 1);
ab   = zeros(68, 1);
ab_p = zeros(68, 1);
for i = 1:68
    M   = EI(:, i);
    med = bootstrap_mediation(zscore(X), zscore(M), zscore(Y), Z_z);
    b(i)    = med.bhat;
    b_p(i)  = med.p_b;
    ab(i)   = med.indirect_hat;
    ab_p(i) = med.p_ab;
end
output_struct.regional_b  = b;
output_struct.regional_ab = ab;

%% Save
if ~exist('../output', 'dir')
    mkdir('../output')
end
save('../output/MACC.mat', 'output_struct')

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
A_design = [Xcol, Zmat];
mdl  = fitglm(A_design, Mcol, 'Distribution', 'normal', 'Link', 'identity');
ahat = mdl.Coefficients.Estimate(2);

% Path b and direct effect c': Y ~ M + X + Z
B_design    = [Mcol, Xcol, Zmat];
mdl         = fitglm(B_design, Ycol, 'Distribution', 'normal', 'Link', 'identity');
bhat        = mdl.Coefficients.Estimate(2);
c_prime_hat = mdl.Coefficients.Estimate(3);

% Total effect c: Y ~ X + Z
C_design = [ones(n,1), Xcol, Zmat];
coefC = regress(Ycol, C_design);
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
    A_design_b = [ones(n,1), Xcol(idx), Zmat(idx,:)];
    B_design_b = [ones(n,1), Mcol(idx), Xcol(idx), Zmat(idx,:)];
    C_design_b = [ones(n,1), Xcol(idx), Zmat(idx,:)];

    coefA_b = regress(Mcol(idx), A_design_b);
    coefB_b = regress(Ycol(idx), B_design_b);
    coefC_b = regress(Ycol(idx), C_design_b);

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
