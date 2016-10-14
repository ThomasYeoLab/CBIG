function [betas, stats] = CBIG_fitGLM_3t(subInfoFile, GMICVFile, q)

% [betas, stats] = CBIG_fitGLM_3t(subInfoFile, GMICVFile, q)
% What is the subtype order?
% Offset by 1, the RID column
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

SUB_COL = 1+3;
COR_COL = 1+1;
fprintf('\n\n\n');
fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');

% All these info is baseline info
[RID_dx, RID_ICVoverGM, RID_prob, RID_age, RID_gender, RID_edu] = CBIG_get_data(subInfoFile, GMICVFile, 3);

% RID of 188 AD, 147 a+ MCI, 43 a+ AD
RID_188AD = RID_dx(RID_dx(:, 2)==3, 1);
RID_147MCI = CBIG_select_predementedGrp(RID_dx, 'a+_mci_147');
RID_43HC = CBIG_select_predementedGrp(RID_dx, 'a+_hc_43');
RID = sort([RID_188AD; RID_147MCI; RID_43HC]);

[qoi, ~] = get_quantityOfInterest(RID, q, 0);
RID_time_score = CBIG_compute_time(qoi, qoi);
RID_time_score = RID_time_score(RID_time_score(:, 2)==0, :); % Baseline only

%--------------------------------------------- GLM

% Construct X
X = zeros(size(RID_time_score, 1), 12);
for idx = 1:size(X, 1)
    RID = RID_time_score(idx, 1);
    % Go fetch this subject's data
    rowIdx = RID_prob(:, 1)==RID;
    dx = RID_dx(rowIdx, 2);
    m = dx==2;
    d = dx==3;
    s = RID_prob(rowIdx, SUB_COL);
    c = RID_prob(rowIdx, COR_COL);
    a = RID_age(rowIdx, 2);
    g = RID_gender(rowIdx, 2);
    e = RID_edu(rowIdx, 2);
    i = RID_ICVoverGM(rowIdx, 2);
    % Fill in one row
    X(idx, :) = [m d s c m*s m*c d*s d*c a g e i];
end

% Construct y
y = RID_time_score(:, 3);

% Fit
[betas, ~, stats] = glmfit(X, y);
