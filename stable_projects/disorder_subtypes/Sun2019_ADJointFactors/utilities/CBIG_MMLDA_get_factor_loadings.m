function rid_prob_reorder = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order)
% rid_prob_reorder = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, order)
%
% Get reordered factor loadings based on rid file, inferred gamma file and the order of factors.
%
% Input:
%   - rid_file      : rid file name
%   - inf_gamma_file: inferred gamma file name
%   - order         : 1 x K array. K is number of factors. e.g. [3, 2, 1] means the 3rd column of gamma file correspond 
%                     to MTL-Memory factor; the 2nd column correspond to Lateral Temporal-Language factor; the 1st column
%                     corresponds to Posterior Cortical-Executive factor in our paper.
%
% Output:
%   - rid_prob_reorder  : N x (K+1) array. N is number of subjects. 1st is rid. 2nd to K+1 th are factor loadings.
%
% Example:
%   rid_prob_reorder = CBIG_MMLDA_get_factor_loadings(rid_file, inf_gamma_file, [3 2 1])
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step3_analyses_internalUse/functions'])

rid = load(rid_file);
prob = CBIG_get_prob(inf_gamma_file);
prob_reorder = prob(:, order);
rid_prob_reorder = [rid prob_reorder];

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Zhang2016_ADFactors/step3_analyses_internalUse/functions'])
