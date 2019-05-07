function edu = CBIG_MMLDA_get_edu(PTDEMOG_file, phase, rid, viscode)
% edu = CBIG_MMLDA_get_edu(PTDEMOG_file, phase, rid, viscode)
%
% Get education year from PTDEMOG spreadsheet.
%
% Input:
%   - PTDEMOG_file  : Demographic spreadsheet 
%   - phase         : cell array of phase for subjects. Each cell contains 'ADNI1', 'ADNIGO' or ADNI2'.
%   - rid           : cell array of rid for subjects
%   - viscode       : For 'ADNI1' subjects, it is cell array of VISCODE, corrsponding to rid.
%                     For 'ADNIGO' or 'ADNI2' subjects, it is cell array of VISCODE2, corresponding to rid. 
%
% Output:
%   - edu           : year of education
%
% Example:
%   PTDMOG_file = [ADNI_spreadsheet_path '/Subject_Characteristics/PTDEMOG.csv'];
%   phase = {'ADNI2', 'ADNI2'};
%   rid = {'3', '23'};
%   viscode = {'sc', 'sc'};
%   edu = CBIG_MMLDA_get_edu(PTDEMOG_file, phase, rid, viscode)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% get edu
fieldnames = {'PTEDUCAT'};
edu = CBIG_MMLDA_get_num_entry(PTDEMOG_file, phase, rid, viscode, fieldnames, 0);