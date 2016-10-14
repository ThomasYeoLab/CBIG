function [betas, stats] = CBIG_fitGLM(subInfoFile, GMICVFile, q)

% [betas, stats] = CBIG_fitGLM(subInfoFile, GMICVFile, q)
% All these info is baseline info
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[RID_dx, RID_ICVoverGM, ~, RID_age, RID_gender, RID_edu] = get_data(subInfoFile, GMICVFile, 3);

% RID of 188 AD, 147 a+ MCI, 43 a+ AD
RID_188AD = RID_dx(RID_dx(:, 2)==3, 1);
RID_147MCI = CBIG_select_predementedGrp(RID_dx, 'a+_mci_147');
RID_43HC = CBIG_select_predementedGrp(RID_dx, 'a+_hc_43');
RID = sort([RID_188AD; RID_147MCI; RID_43HC]);

[qoi, ~] = get_quantityOfInterest(RID, q, 0);
RID_time_score = CBIG_compute_time(qoi, qoi);
RID_time_score = RID_time_score(RID_time_score(:, 2)==0, :); % Baseline only

%--------------------------------------------- GLM

scores_AD = [];
scores_MCI = [];
scores_CN = [];

% Construct X
X = zeros(size(RID_time_score, 1), 6);
for idx = 1:size(X, 1)
    RID = RID_time_score(idx, 1);
    % Go fetch this subject's data
    rowIdx = RID_age(:, 1)==RID;
    dx = RID_dx(rowIdx, 2);
    switch dx
        case 1
            scores_CN = [scores_CN; RID_time_score(idx, 3)];
        case 2
            scores_MCI = [scores_MCI; RID_time_score(idx, 3)];
        otherwise
            scores_AD = [scores_AD; RID_time_score(idx, 3)];
    end
    m = dx==2;
    d = dx==3;
    a = RID_age(rowIdx, 2);
    g = RID_gender(rowIdx, 2);
    e = RID_edu(rowIdx, 2);
    i = RID_ICVoverGM(rowIdx, 2);
    % Fill in one row
    X(idx, :) = [m d a g e i];
end

% Construct y
y = RID_time_score(:, 3);

% Fit
[betas, ~, stats] = glmfit(X, y);

%--------------------------------------------- Report means

fprintf('\n\n');
fprintf('CN mean: %f\n', mean(scores_CN));
fprintf('A+ MCI mean: %f\n', mean(scores_MCI));
fprintf('A+ AD mean: %f\n', mean(scores_AD));