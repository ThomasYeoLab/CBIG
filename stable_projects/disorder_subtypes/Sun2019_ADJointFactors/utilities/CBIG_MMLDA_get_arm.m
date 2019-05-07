function arm = CBIG_MMLDA_get_arm(ARM_file, rid)
% arm = CBIG_MMLDA_get_arm(ARM_file, rid)
%
% Get ARM value and combine it with diagnosis to determine whether it's LMCI, 
% EMCI or SMC.
%
% ARM:
% 1=NL (ADNI1 1.5T only)
% 2=LMCI (ADNI1 1.5T only)
% 3=AD (ADNI1 1.5T only)
% 4=NL (ADNI1 PET+1.5T)
% 5=LMCI (ADNI1 PET+1.5T)
% 6=AD (ADNI1 PET+1.5T)
% 7=NL (ADNI1 3T+1.5T)
% 8=LMCI (ADNI1 3T+1.5T)
% 9=AD (ADNI1 3T+1.5T)
% 10=EMCI
% 11=SMC
%
% 5 Diagnosis:
% 1=NL
% 2=SMC
% 3=EMCI
% 4=LMCI
% 5=AD
%
% mapping between ARM and Diagnosis
% CN and arm!=11 -> NL
% CN and arm==11 -> SMC
% MCI and arm==10 -> EMCI
% MCI and arm!=10 -> LMCI
% AD -> AD
%
% Input:
%   - ARM_file      : path of ARM.csv spreadsheet
%   - rid           : cell array of rid for subjects
%
% Output:
%   - arm           : N x 1 array. N is number of subjects. Each row is number 
%                     ranging from 1 to 11 and the meaning has been show above.
%
% Example:
%   arm = CBIG_MMLDA_get_arm(ARM_file, rid)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

phase = repmat({'ADNIGO2'}, length(rid), 1);
viscode = repmat({'sc'}, length(rid), 1);
fieldname = {'ARM'};
arm = CBIG_MMLDA_get_num_entry(ARM_file, phase, rid, viscode, fieldname, 0, 0);


