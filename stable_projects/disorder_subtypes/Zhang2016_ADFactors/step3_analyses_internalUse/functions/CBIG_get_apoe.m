function apoe = CBIG_get_apoe(rid_of_interest)

% apoe = CBIG_get_apoe(rid_of_interest)
% Fetch APOE from .csv files
%
% Input:
%    - rid_of_interest: the RID's of subjects of interest
%
% Output:
%    - apoe: 1st col is RID, 2nd is VISCODE, 3rd is APOE2, 4th is APOE4

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/APOERES.csv', 'r');
data = textscan(file_id, ...
    '%*q %*q %q %*q %q %*q %*q %q %q %q %*[^\n]', ...
    'delimiter', {',', '\n'}, 'HeaderLines', 1); % read in RID, PTGENDER, and PTDOBYY
fclose(file_id);
data = [data{1, 1} data{1, 2} data{1, 3} data{1, 4} data{1, 5}];

% ADNI-1 only
data = data(1:1159, :);

% Certain RID only
rid_col = cellfun(@str2num, data(:, 1));
keep_ind = ismember(rid_col, rid_of_interest);
data_rid = data(keep_ind, :);

% Sort by RID
data_rid_sort = data_rid;
rid_col = cellfun(@str2num, data_rid(:, 1));
viscode_col = data_rid(:, 2);
examdate_col = data_rid(:, 3);
gen1_col = cellfun(@str2num, data_rid(:, 4));
gen2_col = cellfun(@str2num, data_rid(:, 5));
[rid_col_sort, sort_ind] = sort(rid_col);
data_rid_sort(:, 1) = num2cell(rid_col_sort);
data_rid_sort(:, 2) = viscode_col(sort_ind);
data_rid_sort(:, 3) = examdate_col(sort_ind);
gen1_col_sort = gen1_col(sort_ind);
data_rid_sort(:, 4) = num2cell(gen1_col_sort);
gen2_col_sort = gen2_col(sort_ind);
data_rid_sort(:, 5) = num2cell(gen2_col_sort);

% Overwrite for simplicity
data = data_rid_sort;

rid_gen1_gen2 = cell2mat(data(:, [1 4 5]));
no4 = sum(rid_gen1_gen2(:, [2 3])==4, 2); % How many 4's per row (subject)?
no2 = sum(rid_gen1_gen2(:, [2 3])==2, 2); % How many 2's per row (subject)?

% Add dummy VISCODE columns
apoe = cell(size(data, 1), 4);
apoe(:, 1) = num2cell(rid_gen1_gen2(:, 1));
apoe(:, 2) = repmat({'bl'}, size(no4, 1), 1);
apoe(:, 3) = num2cell(no2);
apoe(:, 4) = num2cell(no4);
