function Z_flip = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, behavior_sign, refpara)
% Z_flip = CBIG_MMLDA_behavior_to_zscore(behavior, regressor, refpara)
%
% Convert behavior matrix to zscore with respect to reference cohort.
% 1. Regress out regressors with respect to reference.
% 2. Zscore withe respect ot reference.
% 3. Flip z score so that larger means worse performance 
%
% Input:
%   - behavior              : N x B matrix, N is # of subjects, B is # of behaviors
%   - regressor             : N x M matrix, N is # of subjects, M is # of regressors
%   - behavior_sign         : cell array with length B. 
%                             If cell value is 'PLUS', bigger means worse performance
%                             If cell value is 'MINUS', smaller means worse performance.
%   - refpara               : output from function CBIG_MMLDA_behavior_to_zscore_create_refpara.m
%
% Output:
%   - Z_flip                : N x B matrix, flipped zscore
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% load reference parameters
load(refpara)

%%%
% regress out nuisance regressors
%%%
% Regress out regressors for all subjects by computing Y - X*b_CN +
% mean(Y_CN) for each voxel
Y = behavior;
X = [ones(size(Y, 1), 1) regressor];
behavior_reg = bsxfun(@plus, Y-X*ref_beta, ref_mean);

%%%
% z score with respect to CN
%%%
Z = bsxfun(@minus, behavior_reg, ref_reg_mean);
Z = bsxfun(@rdivide, Z, ref_reg_std);

%%%
% flip z score
%%%
for i = 1:length(behavior_sign)
    if strcmp(behavior_sign{i}, 'PLUS')
        Z_flip(:, i) = Z(:, i);
    elseif strcmp(behavior_sign{i}, 'MINUS')
        Z_flip(:, i) = -Z(:, i);
    end
end
