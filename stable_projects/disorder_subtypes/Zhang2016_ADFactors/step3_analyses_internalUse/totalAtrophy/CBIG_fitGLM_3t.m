function [betas, stats] = CBIG_fitGLM_3t(subInfoFile, GMICVFile)

% [betas, stats] = CBIG_fitGLM_3t(subInfoFile, GMICVFile)
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

RID_AD = CBIG_select_predementedGrp(RID_dx, 'ad_188');
RID_aMCI = CBIG_select_predementedGrp(RID_dx, 'a+_mci_147');
RID_aHC = CBIG_select_predementedGrp(RID_dx, 'a+_hc_43');
RID = [RID_AD; RID_aMCI; RID_aHC];

%--------------------------------------------- GLM

% Construct X and y
X = zeros(numel(RID), 5);
y = zeros(numel(RID), 1);
for idx = 1:size(X, 1)
    % Go fetch this subject's data
    rowIdx = RID_prob(:, 1)==RID(idx);
    s = RID_prob(rowIdx, SUB_COL);
    c = RID_prob(rowIdx, COR_COL);
    a = RID_age(rowIdx, 2);
    g = RID_gender(rowIdx, 2);
    e = RID_edu(rowIdx, 2);
    % Fill in one row
    X(idx, :) = [s c a g e];
    y(idx) = RID_ICVoverGM(rowIdx, 2);
end

% Fit
[betas, ~, stats] = glmfit(X, y);
