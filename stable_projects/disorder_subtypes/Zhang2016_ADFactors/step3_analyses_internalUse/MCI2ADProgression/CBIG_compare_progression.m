function CBIG_compare_progression(subInfoFile, GMICVFile, windowSize, k)

% CBIG_compare_progression(subInfoFile, GMICVFile, windowSize, k)
% What is the subtype order?
% Offset by 1, the RID column
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

TEM_COL = 1+2;
SUB_COL = 1+3;
fprintf('\n\n\n');
fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');

%------------------------------- Convert or not

[rid_diagnosis, rid_ICVoverGM, rid_prob, rid_age, rid_gender, rid_edu] = CBIG_get_data(subInfoFile, GMICVFile, k);
% ridOfInterest = CBIG_select_predementedGrp(rid_diagnosis, 'mci_394');
ridOfInterest = CBIG_select_predementedGrp(rid_diagnosis, 'a+_hcridOfInterest = select_predementedGrp(rid_diagnosis, 'a+_hc&mci_190');mci_190');
% ridOfInterest = CBIG_select_predementedGrp(rid_diagnosis, 'a+_mci_147');

[rid_conv, rid_nonconv] = CBIG_find_conv_nonconv(ridOfInterest, windowSize);

fprintf('# convertors = %i\n', numel(rid_conv));
fprintf('# nonconvertors = %i\n', numel(rid_nonconv));
fprintf('# na = %i\n', numel(ridOfInterest)-numel(rid_conv)-numel(rid_nonconv));
fprintf('# convertors + nonconvertors = %i\n', numel(rid_conv)+numel(rid_nonconv));
fprintf('\n');

%!!!!!!!!! 0 for nonconv, 1 for conv
% x-axis is flipped
% such that left is always "former is bad"
rid_convValue = repmat(ridOfInterest, [1, 2]);
conv_ind = ismember(rid_convValue(:, 1), rid_conv);
rid_convValue(conv_ind, 2) = 1; % conv
nonconv_ind = ismember(rid_convValue(:, 1), rid_nonconv);
rid_convValue(nonconv_ind, 2) = 0; % nonconv
rid_convValue(~(conv_ind|nonconv_ind), :) = []; % neither

%------------------------------- Logistic regression

X = zeros(size(rid_convValue, 1), 6);
for idx = 1:size(X, 1)
    rid = rid_convValue(idx, 1);
    temProb = rid_prob(rid_prob(:, 1)==rid, TEM_COL);
    subProb = rid_prob(rid_prob(:, 1)==rid, SUB_COL);
    age = rid_age(rid_age(:, 1)==rid, 2);
    gender = rid_gender(rid_gender(:, 1)==rid, 2);
    edu = rid_edu(rid_edu(:, 1)==rid, 2);
    ICVoverGM = rid_ICVoverGM(rid_ICVoverGM(:, 1)==rid, 2);
    X(idx, :) = [temProb subProb age gender edu ICVoverGM];
end
y = rid_convValue(:, 2);
n = ones(size(y, 1), 1); % each subject, one trial

[b, ~, stats] = glmfit(X, [y n], 'binomial', 'link', 'logit'); % log(�/(1-�)) = Xb

%------------------------------- Compare likelihood

% Compute the conversion likelihood for "average" male and female
avg_stats = [mean(X(:, 3)) -1 mean(X(:, 5)) mean(mean(X(:, 6)))];
convLike_tem_m = CBIG_compute_convLike(b, [1 0 0], 1, avg_stats);
convLike_tem_f = CBIG_compute_convLike(b, [1 0 0], 2, avg_stats);
convLike_sub_m = CBIG_compute_convLike(b, [0 1 0], 1, avg_stats);
convLike_sub_f = CBIG_compute_convLike(b, [0 1 0], 2, avg_stats);
convLike_cor_m = CBIG_compute_convLike(b, [0 0 1], 1, avg_stats);
convLike_cor_f = CBIG_compute_convLike(b, [0 0 1], 2, avg_stats);

% How many times more likely?
fprintf('F - temporal/subcortical = %f\n', convLike_tem_f/convLike_sub_f);
fprintf('M - temporal/subcortical = %f\n', convLike_tem_m/convLike_sub_m);
fprintf('F - temporal/cortical = %f\n', convLike_tem_f/convLike_cor_f);
fprintf('M - temporal/cortical = %f\n', convLike_tem_m/convLike_cor_m);
fprintf('\n');

%--------------------------------- Hypothesis testing: exact LR test

disp('Exact LR test');

% Unrestricted model
b_u = b;
yFit = glmval(b_u, X, 'logit', 'size', n);
ll_u = sum(log(binopdf(y, n, yFit./n)));

A = zeros(3, 2);
% Forest plot data
% S - T
H = [0 -1 1 0 0 0 0];
idx1 = find(H==1);
idx2 = find(H==-1);
[mu, sigma2] = CBIG_compute_mu_sigma2(stats, idx1, idx2);
A(1, :) = [mu, sqrt(sigma2)];
% C - T
H = [0 -1 0 0 0 0 0];
idx1 = find(H==-1);
A(2, :) = [-stats.beta(idx1), sqrt(stats.covb(idx1, idx1))];
% C - S
H = [0 0 -1 0 0 0 0];
idx1 = find(H==-1);
A(3, :) = [-stats.beta(idx1), sqrt(stats.covb(idx1, idx1))];

% Overall
% Restricted model
X_r = [X(:, 3) X(:, 4) X(:, 5) X(:, 6)];
dof = 2;
disp('Overall');
p = CBIG_lr_test(X_r, y, n, ll_u, dof)

% S - T
% Restricted model
X_r = [X(:, 1)+X(:, 2) X(:, 3) X(:, 4) X(:, 5) X(:, 6)];
dof = 1;
disp('S - T');
p = CBIG_lr_test(X_r, y, n, ll_u, dof)

% C - T
% Restricted model
X_r = [X(:, 2) X(:, 3) X(:, 4) X(:, 5) X(:, 6)];
dof = 1;
disp('C - T');
p = CBIG_lr_test(X_r, y, n, ll_u, dof)

% C - S
% Restricted model
X_r = [X(:, 1) X(:, 3) X(:, 4) X(:, 5) X(:, 6)];
dof = 1;
disp('C - S');
p = CBIG_lr_test(X_r, y, n, ll_u, dof)

% Forest plot
close all;
CBIG_forestPlot_xRev(A);

% %--------------------------------- Hypothesis testing: approximate linear hypothesis test
% 
% %---------------- Overall
% 
% H = [0 1 -1 0 0 0 0; 0 1 0 0 0 0 0];
% c = [0; 0];
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
% 
% %---------------- Pairwise comparisons
% 
% A = zeros(3, 2);
% 
% % T vs. S
% H = [0 1 -1 0 0 0 0];
% c = 0;
% disp('T vs. S');
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
% idx1 = find(H==1);
% idx2 = find(H==-1);
% [mu, sigma2] = CBIG_compute_mu_sigma2(stats, idx1, idx2);
% A(1, :) = [mu, sqrt(sigma2)];
% 
% % T vs. C
% H = [0 1 0 0 0 0 0];
% c = 0;
% disp('T vs. C');
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
% idx1 = find(H==1);
% A(2, :) = [stats.beta(idx1), sqrt(stats.covb(idx1, idx1))];
% 
% % S vs. C
% H = [0 0 1 0 0 0 0];
% c = 0;
% disp('S vs. C');
% p = linhyptest(b, stats.covb, c, H, stats.dfe)
% idx1 = find(H==1);
% A(3, :) = [stats.beta(idx1), sqrt(stats.covb(idx1, idx1))];
