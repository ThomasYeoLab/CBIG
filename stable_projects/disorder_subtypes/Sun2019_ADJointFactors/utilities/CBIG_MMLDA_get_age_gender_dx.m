function [age, gender, dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase, rid, viscode)
% [age, gender, dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file,...
%   phase, rid, viscode)
%
% This function get age, gender, dx based on ADNI DXSUM and PTDEMOG spreadsheets. 
%
% Input:
%   - DXSUM_file    : Diagnosis summary spreadsheet
%   - PTDEMOG_file  : Demographic spreadsheet 
%   - phase         : cell array of phase for subjects. Each cell contains 'ADNI1', 'ADNIGO' or ADNI2'.
%   - rid           : cell array of rid for subjects
%   - viscode       : For 'ADNI1' subjects, it is cell array of VISCODE, corrsponding to rid.
%                     For 'ADNIGO' or 'ADNI2' subjects, it is cell array of VISCODE2, corresponding to rid. 
% 
% Output:
%   - age           : a N x 1 array. N is the # of subjects.
%   - gender        : a N x 1 array. N is the # of subjects.
%   - dx            : a N x 1 array. 1 means CN, 2 means MCI, 3 means AD.
%
% Example:
% ADNI_speadsheet_path = '/share/users/imganalysis/yeolab/data/ADNI/All/ADNI_161017/documentation';
% DXSUM_file = [ADNI_spreadsheet_path '/Assessments/DXSUM_PDXCONV_ADNI_ALL.csv'];
% PTDMOG_file = [ADNI_spreadsheet_path '/Subject_Characteristics/PTDEMOG.csv'];
% phase = {'ADNI2', 'ADNI2'};
% rid = {'3', '23'};
% viscode = {'bl', 'm06'};
% [age, gender, dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase, rid, viscode)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% get dx
fieldnames = {'DXCURREN', 'DXCHANGE', 'DIAGNOSIS'};
dxcurr_dxchange_diagnosis = CBIG_MMLDA_get_num_entry(DXSUM_file, phase, rid, viscode, fieldnames, 0);
dx = zeros(length(dxcurr_dxchange_diagnosis), 1);
for i = 1:length(dxcurr_dxchange_diagnosis)
    if ~isnan(dxcurr_dxchange_diagnosis(i, 1))
        dx(i) = dxcurr_dxchange_diagnosis(i, 1);
    elseif ~isnan(dxcurr_dxchange_diagnosis(i, 2))
        dx(i) = CBIG_MMLDA_dxchange2dxcurr(dxcurr_dxchange_diagnosis(i, 2));
    else
        dx(i) = dxcurr_dxchange_diagnosis(i, 3);
    end
end

% get examdate
fieldnames = {'EXAMDATE'};
examdate = CBIG_MMLDA_get_string_entry(DXSUM_file, phase, rid, viscode, fieldnames, 0);
examyy_exammm = zeros(length(examdate), 2);
for i = 1:length(examdate)
    yy_mm_dd = strsplit(examdate{i}, '-');
    if isnan(str2num(yy_mm_dd{1}))
        examyy_exammm(i, 1) = NaN;
        examyy_exammm(i, 2) = NaN;
    else
        examyy_exammm(i, 1) = str2num(yy_mm_dd{1});
        examyy_exammm(i, 2) = str2num(yy_mm_dd{2});
    end
end

% get gender, year of birth, month of birth
fieldnames = {'PTGENDER', 'PTDOBYY', 'PTDOBMM'};
gender_dobyy_dobmm = CBIG_MMLDA_get_num_entry(PTDEMOG_file, phase, rid, repmat({'sc'}, length(rid), 1), fieldnames, 0, 0);
gender = gender_dobyy_dobmm(:, 1);

% calculate the age
age = examyy_exammm(:, 1) - gender_dobyy_dobmm(:, 2) + (examyy_exammm(:, 2) - gender_dobyy_dobmm(:, 3)) / 12;


