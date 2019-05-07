function [score, state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, ...
    UCBERKELEYAV45_file, amyloid_phase, phase, rid, viscode)
% [score, state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, amyloid_phase, phase, rid, viscode)
%
% Get amyloid score and state from UPENNBIOMK and UCBERKELEYAV45 spreadsheets.
%
% Input:
%   - UPENNBIOMK_file       : path of UPENNBIOMK.csv spreadsheet 
%   - UCBERKELEYAV45_file   : path of UCBERKELEYAV45.csv spreadsheet
%   - amyloid_phase         : if 'ADNI1', you use CSF amyloid from UPENNBIOMK
%                             if 'ADNI2', you use PET amyloid from UCBERKELEYAV45
%   - phase                 : cell array of phase for subjects. Each cell contains 'ADNI1', 'ADNIGO' or ADNI2'.
%   - rid                   : cell array of rid for subjects
%   - viscode               : For 'ADNI1' subjects, it is cell array of VISCODE, corrsponding to rid.
%                             For 'ADNIGO' or 'ADNI2' subjects, it is cell array of VISCODE2, corresponding to rid. 
%
% Output:
%   - score                 : N x 1 array. N is number of subjects. CSF pg/ml or PET SUVR.
%   - state                 : N x 1 array. 1 means amyloid beta positive. 0 means amyloid beta negative.
%
% Example:
%   [score, state] = CBIG_MMLDA_get_amyloid(UPENNBIOMK_file, UCBERKELEYAV45_file, amyloid_phase, phase, rid, viscode)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if strcmp(amyloid_phase, 'ADNI1')
    % read CSF Amyloid from spreadsheet
    score = CBIG_MMLDA_get_num_entry(UPENNBIOMK_file, phase, rid, viscode, {'ABETA142'}, 0);
    state = nan(length(score), 1);
    state(score < 192) = 1;
    state(score >= 192) = 0;
elseif strcmp(amyloid_phase, 'ADNI2')
    % read PET Amyloid from spreadsheet 
    score = CBIG_MMLDA_get_num_entry(UCBERKELEYAV45_file, phase, rid, viscode, ...
        {'SUMMARYSUVR_WHOLECEREBNORM'}, 0);
    state = CBIG_MMLDA_get_num_entry(UCBERKELEYAV45_file, phase, rid, viscode, ...
        {'SUMMARYSUVR_WHOLECEREBNORM_1_11CUTOFF'}, 0);
end