function  T_S_C_p = CBIG_gender(subInfoFile, GMICVFile, K)

% T_S_C_p = CBIG_gender(subInfoFile, GMICVFile, K)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

T_S_C_p = cell(1, 4);


display('************* SEX *************');
fprintf('\n');

%------------------------------- Subtypes

% Subtypes estimated at ADNI 1 baseline
[rid_dx, ~, rid_prob, ~, rid_gender, ~] = get_data(subInfoFile, GMICVFile, K);
rid_prob = rid_prob(rid_dx(:, 2)==3, :);
rid_gender = rid_gender(rid_dx(:, 2)==3, :);

%------------------------------- Compute stats

% 1-1=0: male; 2-1=1: female
y = rid_gender(:, 2)-1;
[wMean, wStdDev] = CBIG_compute_weightedMeanStdDev(rid_prob(:, 2:end), y);

if K == 3
    % What is the subtype order?
    % Offset by 1, the RID column
    TEM_COL = 1+2;
    SUB_COL = 1+3;
    COR_COL = 1+1;
    fprintf('\n\n\n');
    fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
    fprintf('!!! Make sure this is true in the gamma file.\n');
    %--------------------------------- Hypothesis testing: exact LR test
    X = rid_prob(:, [TEM_COL SUB_COL]);
    n = ones(size(y, 1), 1);
    [b, ~, stats] = glmfit(X, [y n], 'binomial', 'link', 'logit'); % log(�/(1-�)) = Xb
    disp('Exact LR test');
    % Unrestricted model
    b_u = b;
    yFit = glmval(b_u, X, 'logit', 'size', n);
    ll_u = sum(log(binopdf(y, n, yFit./n)));
    % Overall
    % Restricted model
    X_r = [];
    dof = 2;
    disp('Overall');
    p = CBIG_lr_test(X_r, y, n, ll_u, dof)
    % Restricted model
    X_r = [X(:, 1)+X(:, 2)];
    dof = 1;
    disp('T == S?');
    CBIG_lr_test(X_r, y, n, ll_u, dof)
    % Restricted model
    X_r = [X(:, 2)];
    dof = 1;
    disp('T == C?');
    CBIG_lr_test(X_r, y, n, ll_u, dof)
    % Restricted model
    X_r = [X(:, 1)];
    dof = 1;
    disp('S == C?');
    CBIG_lr_test(X_r, y, n, ll_u, dof)
else
    error('Unconfigured; refer to older versions.');
end

T_S_C_p{1} = sprintf('%.3f (%.3f)', wMean(TEM_COL-1), wStdDev(TEM_COL-1));
T_S_C_p{2} = sprintf('%.3f (%.3f)', wMean(SUB_COL-1), wStdDev(SUB_COL-1));
T_S_C_p{3} = sprintf('%.3f (%.3f)', wMean(COR_COL-1), wStdDev(COR_COL-1));
T_S_C_p{4} = sprintf('%e', p);

% %--------------------------------- Hypothesis testing: approximate linear hypothesis test
%
% %---------------- Overall
%
% disp('Overall');
% % For linear regression, the p-values are exact.
% % H*b = c
% H = [0 1 -1; 0 1 0];
% c = [0; 0];
% % The dfd is often called the degrees of freedom error or dfe. In the simplest case of a one-factor between-subjects ANOVA,
% % dfn = a-1
% % dfd = N-a
% % where "a" is the number of groups and "N" is the total number of subjects in the experiment.
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
%
% %---------------- Pairwise comparisons
%
% % T vs. S
% H = [0 1 -1];
% c = 0;
% disp('T vs. S');
% if H*b > 0
%     fprintf('T > S\n');
% else
%     fprintf('T < S\n');
% end
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
%
% % T vs. C
% H = [0 1 0];
% c = 0;
% disp('T vs. C');
% if H*b > 0
%     fprintf('T > C\n');
% else
%     fprintf('T < C\n');
% end
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
%
% % S vs. C
% H = [0 0 1];
% c = 0;
% disp('S vs. C');
% if H*b > 0
%     fprintf('S > C\n');
% else
%     fprintf('S < C\n');
% end
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
