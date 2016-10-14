function ps = CBIG_pairwiseCmpMeans(subInfoFile, GMICVFile, Q)

% ps = CBIG_pairwiseCmpMeans(subInfoFile, GMICVFile, Q)
% All these info is baseline info
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[RID_dx, ~, ~, ~, ~, ~] = get_data(subInfoFile, GMICVFile, 3);

% RIDs of 188 AD, 147 a+ MCI, 43 a+ CN
RID_188AD = RID_dx(RID_dx(:, 2)==3, 1);
RID_147MCI = CBIG_select_predementedGrp(RID_dx, 'a+_mci_147');
RID_43CN = CBIG_select_predementedGrp(RID_dx, 'a+_hc_43');
RID = sort([RID_188AD; RID_147MCI; RID_43CN]);

% Fetch scores
[qoi, ~] = get_quantityOfInterest(RID, Q, 0);
RID_time_score = CBIG_compute_time(qoi, qoi);
RID_score = RID_time_score(RID_time_score(:, 2)==0, [1 3]); % Baseline only

%% One-way ANOVA
n = size(RID_score, 1);
y = zeros(1, n);
grp = cell(1, n);
for idx = 1:n
    % y
    y(idx) = RID_score(idx, 2);
    % grp
    pt = RID_score(idx, 1);
    if sum(ismember(RID_188AD, pt)) == 1
        grp{idx} = 'AD';
    elseif sum(ismember(RID_147MCI, pt)) == 1
        grp{idx} = 'A+ MCI';
    elseif sum(ismember(RID_43CN, pt)) == 1
        grp{idx} = 'A+ CN';
    else
        error('Wrong!');
    end
end
[p, tbl, stats] = anova1(y, grp);

% % Multiple comparisons based on stats from ANOVA
% [c, m, ~, gnames] = multcompare(stats);
% [{'GrpName' 'Mean' 'SE'};...
%     gnames num2cell(m)]
% [{'Grp1' 'Grp2' '95% LowLim' 'MeanDiff' '95% UppLim' 'p'};...
%     gnames(c(:, 1)) gnames(c(:, 2)) num2cell(c(:, 3:6))]

%% Two-sample t-tests

% AD
ind = ismember(RID_score(:, 1), RID_188AD);
scores_AD = RID_score(ind, 2);

% MCI
ind = ismember(RID_score(:, 1), RID_147MCI);
scores_MCI = RID_score(ind, 2);

% CN
ind = ismember(RID_score(:, 1), RID_43CN);
scores_CN = RID_score(ind, 2);

% Means
fprintf('\nAD mean = %f\n', mean(scores_AD));
fprintf('MCI mean = %f\n', mean(scores_MCI));
fprintf('CN mean = %f\n', mean(scores_CN));

ps = [];

% AD vs. MCI
[~, p] = ttest2(scores_AD, scores_MCI);
fprintf('\nAD vs. MCI: p = %e\n', p);
ps = [ps; p];

% AD vs. CN
[~, p] = ttest2(scores_AD, scores_CN);
fprintf('AD vs. CN: p = %e\n', p);
ps = [ps; p];

% MCI vs. CN
[~, p] = ttest2(scores_MCI, scores_CN);
fprintf('MCI vs. CN: p = %e\n', p);
ps = [ps; p];