function mem_ef = CBIG_get_mem_ef(rid, predementedOnly)

% mem_ef = CBIG_get_mem_ef(rid, predementedOnly)
% Fetch MEM and EF scores from .csv files
%
% Input:
%    - rid: the RID's of subjects of interest
%
% Output:
%    - mem_ef: 1st col is RID, 2nd is VISCODE, 3rd is EXAMDATE, 4th is MEM,
%    5th is EF

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% % Ask users if they want raw or re-normalized (w.r.t. baseline of 228 CN) Mem/EF scores
% s = input('\n\n!!!!!!\nRe-normalize Mem/EF scores w.r.t. baseline of 228 CN?\nAnswer y or n\n', 's');
warning('Use raw ADNI-Mem and ADNI-EF by default. Want re-normalized scores? Edit here.');
s = 'n';

if strcmp(s, 'n')
    % Read MMSE.csv
    file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/UWNPSYCHSUM_06_09_14.csv', 'r');
    data = textscan(file_id, ...
        '%q %q %*q %q %*q %q %q %*[^\n]', ...
        'delimiter', {',', '\n'}, 'HeaderLines', 1);
    fclose(file_id);
    data = [data{1, 1} data{1, 2} data{1, 3} data{1, 4} data{1, 5}];
    
    % Remove rows who have an empty MEM or EF score
    data(sum(cellfun(@isempty, data), 2)~=0, :) = [];
    
    % Certain RID only
    rid_col = cellfun(@str2num, data(:, 1));
    keep_ind = ismember(rid_col, rid);
    data_rid = data(keep_ind, :);
    
    % Sort by RID
    data_rid_sort = data_rid;
    rid_col = cellfun(@str2num, data_rid(:, 1));
    viscode_col = data_rid(:, 2);
    examdate_col = data_rid(:, 3);
    mem_col = cellfun(@str2num, data_rid(:, 4));
    ef_col = cellfun(@str2num, data_rid(:, 5));
    [rid_col_sort, sort_ind] = sort(rid_col);
    mem_col_sort = mem_col(sort_ind);
    ef_col_sort = ef_col(sort_ind);
    data_rid_sort(:, 1) = num2cell(rid_col_sort);
    data_rid_sort(:, 2) = viscode_col(sort_ind);
    data_rid_sort(:, 3) = examdate_col(sort_ind);
    data_rid_sort(:, 4) = num2cell(mem_col_sort);
    data_rid_sort(:, 5) = num2cell(ef_col_sort);
    
    mem_ef = data_rid_sort;
    
    if predementedOnly == 0;
        disp('Use all the data we have');
    else
        error('CBIG_get_mem_ef() not fully configured with the "predementedOnly" flag yet');
    end
elseif strcmp(s, 'y')
    load('../normalizeMemEF/normalizedMemEF.mat');
    ridList = cell2mat(rid_viscode_examdate_mem_ef_normalized(:, 1));
    ind = ismember(ridList, rid);
    mem_ef = rid_viscode_examdate_mem_ef_normalized(ind, :);
else
    error('Answer y or n only!');
end