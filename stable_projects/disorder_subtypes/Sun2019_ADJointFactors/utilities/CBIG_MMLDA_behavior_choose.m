function [ADAS, MMSE, NEUROBAT] = CBIG_MMLDA_behavior_choose(phase)
% [ADAS, MMSE, NEUROBAT] = CBIG_MMLDA_behavior_choose(phase)
%
% This function reads in the ADNI1_behavior_choose.csv or ADNI2_behavior_choose.csv
% into matlab and get the filednames, description, sign, category of ADAS, MMSE
% and NEUROBAT behavior scores. For the sign, 'PLUS' means bigger score = worse performance.
% 'MINUS' means smaller score = worse perfermance. For the category, 'MEM' means memory
% 'EF' means executive function, 'NONE' means unknown.
%
% Input:
%   - phase : 'ADNI1' or 'ADNI2' 
%
% Output:
%   - ADAS      : a structre contains FIELDNAMES, DESCRIPTION, SIGN, CATE of ADAS behavior scores
%   - MMSE      : a structre contains FIELDNAMES, DESCRIPTION, SIGN, CATE of MMSE behavior scores
%   - NEUROBAT  : a structre contains FIELDNAMES, DESCRIPTION, SIGN, CATE of NEUROBAT behavior scores
%
% Example:
%   [ADAS, MMSE, NEUROBAT] = CBIG_MMLDA_behavior_choose('ADNI2')
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ADNI_behavior_choose_file = [phase '_behavior_choose.csv'];

% read the behavior choose table
t_choose = readtable(ADNI_behavior_choose_file);

% read the fieldname of each behavioral score
ADAS.FIELDNAMES = t_choose.ADAS(~cellfun(@isempty, t_choose.ADAS));
MMSE.FIELDNAMES = t_choose.MMSE(~cellfun(@isempty, t_choose.MMSE));
NEUROBAT.FIELDNAMES = t_choose.NEUROBAT(~cellfun(@isempty, t_choose.NEUROBAT));

% read the description of each behavioral score
ADAS.DESCRIPTION = t_choose.ADAS_DESCRIPTION(~cellfun(@isempty, t_choose.ADAS_DESCRIPTION));
MMSE.DESCRIPTION = t_choose.MMSE_DESCRIPTION(~cellfun(@isempty, t_choose.MMSE_DESCRIPTION));
NEUROBAT.DESCRIPTION = t_choose.NEUROBAT_DESCRIPTION(~cellfun(@isempty, t_choose.NEUROBAT_DESCRIPTION));

% read the sign of each behavioral score 
% 'PLUS': bigger = worse performance
% 'MINUS': smaller = worse perfermance
ADAS.SIGN = t_choose.ADAS_SIGN(~cellfun(@isempty, t_choose.ADAS_SIGN));
MMSE.SIGN = t_choose.MMSE_SIGN(~cellfun(@isempty, t_choose.MMSE_SIGN));
NEUROBAT.SIGN = t_choose.NEUROBAT_SIGN(~cellfun(@isempty, t_choose.NEUROBAT_SIGN));

% read the category of each behavioral score based on judgement of Thomas and Nanbo
% 'MEM': memory
% 'EF' : executive function
% 'NONE': unknown
ADAS.CATE = t_choose.ADAS_CATE(~cellfun(@isempty, t_choose.ADAS_CATE));
MMSE.CATE = t_choose.MMSE_CATE(~cellfun(@isempty, t_choose.MMSE_CATE));
NEUROBAT.CATE = t_choose.NEUROBAT_CATE(~cellfun(@isempty, t_choose.NEUROBAT_CATE));
