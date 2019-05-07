function behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, NEUROBAT_spreadsheet, ADAS_phase, phase, rid, viscode)
% behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
%     ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, NEUROBAT_spreadsheet, phase, rid, viscode)
%
% This function get ADAS MMSE NEUROBAT scores based on ADNI1_behavior_choose.csv
% or ADNI2_behavior_choose.csv and concat them together.
%
% Input:
%   - ADAS_ADNI1_spreadsheet    : ADNI1 ADAS spreadsheet
%   - ADAS_ADNI2_spreadsheet    : ADNI2 ADAS spreadsheet
%   - MMSE_spreadsheet          : MMSE spreadsheet 
%   - NEUROBAT_spreadsheet      : NEUROBAT spreadsheet
%   - phase         : string 'ADNI1' or ADNI2'.
%   - rid           : cell array of rid for subjects
%   - viscode       : For 'ADNI1' subjects, it is cell array of VISCODE, corrsponding to rid.
%                     For 'ADNIGO' or 'ADNI2' subjects, it is cell array of VISCODE2, corresponding to rid. 
% 
% Output:
%   - behavior_scores           : a N x M matrix. N is the # of subjects and M is # of behaviors.
%
% Example:
% ADNI_spreadsheet_path = '/share/users/imganalysis/yeolab/data/ADNI/All/ADNI_161017/documentation';
% ADAS_ADNI1_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADASSCORES.csv'];
% ADAS_ADNI2_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADAS_ADNIGO2.csv'];
% MMSE_spreadsheet = [ADNI_spreadsheet_path '/Assessments/MMSE.csv'];
% NEUROBAT_spreadsheet = [ADNI_spreadsheet_path '/Assessments/NEUROBAT.csv'];
% phase = 'ADNI2';
% rid = {'3', '23'};
% viscode = {'bl', 'm06'};
% behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
% ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, NEUROBAT_spreadsheet, phase, rid, viscode)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[ADAS_ADNI1, MMSE_ADNI1, NEUROBAT_ADNI1] = CBIG_MMLDA_behavior_choose('ADNI1');
[ADAS_ADNI2, MMSE_ADNI2, NEUROBAT_ADNI2] = CBIG_MMLDA_behavior_choose('ADNI2');

% get behavior scores
ADAS_scores = zeros(length(rid), length(ADAS_ADNI2.FIELDNAMES));
MMSE_scores = zeros(length(rid), length(MMSE_ADNI2.FIELDNAMES));
NEUROBAT_scores = zeros(length(rid), length(NEUROBAT_ADNI2.FIELDNAMES));
if strcmp(ADAS_phase, 'ADNI1')
    ADAS_scores  = CBIG_MMLDA_get_num_entry(ADAS_ADNI1_spreadsheet, repmat({'ADNI1'}, length(rid), 1), rid, viscode, ADAS_ADNI1.FIELDNAMES, 0);
    MMSE_scores = CBIG_MMLDA_get_num_entry(MMSE_spreadsheet, repmat({'ADNI1'}, length(rid), 1), rid, viscode, MMSE_ADNI1.FIELDNAMES, 0);
    NEUROBAT_scores = CBIG_MMLDA_get_num_entry(NEUROBAT_spreadsheet, repmat({'ADNI1'}, length(rid), 1), rid, viscode, NEUROBAT_ADNI1.FIELDNAMES, 0);
elseif strcmp(ADAS_phase, 'ADNI2') || strcmp(ADAS_phase, 'ADNI3')
    ADAS_scores  = CBIG_MMLDA_get_num_entry(ADAS_ADNI2_spreadsheet, phase, rid, viscode, ADAS_ADNI2.FIELDNAMES, 0);
    MMSE_scores = CBIG_MMLDA_get_num_entry(MMSE_spreadsheet, phase, rid, viscode, MMSE_ADNI2.FIELDNAMES, 0);
    NEUROBAT_scores = CBIG_MMLDA_get_num_entry(NEUROBAT_spreadsheet, phase, rid, viscode, NEUROBAT_ADNI2.FIELDNAMES, 0);
end

% set negative value to NaN
ADAS_scores(ADAS_scores < 0)            = NaN;
MMSE_scores(MMSE_scores < 0)            = NaN;
NEUROBAT_scores(NEUROBAT_scores < 0)    = NaN;

% For MMSE we get derivative scores which is based on wiki: MMSE
MMSE_0 = MMSE_scores;
MMSE_0(MMSE_0 == 2) = 0;
MMSE_summary = zeros(size(MMSE_scores, 1), 8);
MMSE_summary(:, 1) = sum(MMSE_0(:, 1:5), 2);
MMSE_summary(:, 2) = sum(MMSE_0(:, 6:10), 2);
MMSE_summary(:, 3) = sum(MMSE_0(:, 11:13), 2);
MMSE_summary(:, 4) = sum(MMSE_0(:, 14:18), 2);
MMSE_summary(:, 5) = sum(MMSE_0(:, 19:21), 2);
MMSE_summary(:, 6) = sum(MMSE_0(:, 22:23), 2);
MMSE_summary(:, 7) = MMSE_0(:, 24);
MMSE_summary(:, 8) = sum(MMSE_0(:, 25:30), 2);

% put ADAS, MMSE, NEUROBAT scores together
behavior_scores = [ADAS_scores MMSE_summary NEUROBAT_scores];