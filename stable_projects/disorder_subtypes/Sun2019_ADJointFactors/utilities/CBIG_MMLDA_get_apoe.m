function [apoe2, apoe4] = CBIG_MMLDA_get_apoe(APOERES_file, rid)
% [apoe2, apoe4] = CBIG_MMLDA_get_apoe(APOERES_file, rid)
%
% Get apoe2 and apoe4 gene number for each subject. 
%
% Input:
%   - APOERES_file  : path of APOERES.csv spreadsheet
%   - rid           : cell array of rid for subjects
%
% Output:
%   - apoe2         : N x 1 array. N is number of subjects. Each row is number of apoe2 genes. e.g. 0, 1, 2
%   - apoe4         : N x 1 array. N is number of subjects. Each row is number of apoe4 genes. e.g. 0, 1, 2
%
% Example:
%   [apoe2, apoe4] = CBIG_MMLDA_get_apoe(APOERES_file, rid)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

phase = repmat({'ADNIGO2'}, length(rid), 1);
viscode = repmat({'sc'}, length(rid), 1);
fieldnames = {'APGEN1', 'APGEN2'};
gen1_gen2 = CBIG_MMLDA_get_num_entry(APOERES_file, phase, rid, viscode, fieldnames, 0, 0);
ind_nonan = ~isnan(sum(gen1_gen2, 2));
apoe4 = nan(length(rid), 1);
apoe2 = nan(length(rid), 1);
apoe4(ind_nonan) = sum(gen1_gen2(ind_nonan, :)==4, 2); % How many 4's per row (subject)?
apoe2(ind_nonan) = sum(gen1_gen2(ind_nonan, :)==2, 2); % How many 2's per row (subject)?