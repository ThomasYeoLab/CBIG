function new_score_names = CBIG_ASDf_renameBehavScores(score_names)
% new_score_names = CBIG_ASDf_renameBehavScores(score_names)
%
% Utils function to rename behavioral scores, mainly for plotting purpose.
% 1. Trim scale name e.g. 'RBSR_6SUBSCALE', 'SRS';
% 2. Trim 'RAW' (raw scores) or 'T' (T scores) characters
% 
% Input:
%     - score_names:
%           M-dim cell array, where each cell contains a string of a
%           behavioral score name, M is the number of scores
%
% Output:
%     - new_score_names:
%           M-dim cell array, renamed behavioral scores
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

score_names = strrep(score_names,'RBSR_6SUBSCALE_','');
score_names = strrep(score_names,'SRS_','');
score_names = strrep(score_names,'_RAW','');
score_names = strrep(score_names,'BRIEF_','');
score_names = strrep(score_names,'CBCL_6_18_','');
new_score_names = strrep(score_names,'_T','');
