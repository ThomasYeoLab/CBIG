function edu = CBIG_get_edu(rid_of_interest)

% edu = CBIG_get_edu(rid_of_interest)
% Fetch education years from .csv files
%
% Input:
%    - rid_of_interest: the RID's of subjects of interest
%
% Output:
%    - edu: 1st col is RID, 2nd col is (dummy) VISCODE, 3rd is (dummy) EXAMDATE, 4th is year

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Read PTDEMOG.csv
file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/PTDEMOG.csv', 'r');
data = textscan(file_id, ...
    '%q %*q %q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %q %*[^\n]',...
    'delimiter', {',', '\n'}, 'HeaderLines', 1);
fclose(file_id);
data = [data{1, 1} data{1, 2} data{1, 3}];

% ADNI-1 only
adni1_ind = ismember(data(:, 1), 'ADNI1');
data = data(adni1_ind, :);

rid_col = cellfun(@str2num, data(:, 2));
keep_ind = ismember(rid_col, rid_of_interest);
data_rid = data(keep_ind, :);

% Remove rows with empty years
keep_ind = ~cellfun(@isempty, data_rid(:, 3));
data_rid = data_rid(keep_ind, :);

% Col 1: RID; Col 2: education
rid = cellfun(@str2num, data_rid(:, 2));
edu = cellfun(@str2num, data_rid(:, 3));
rid_edu = [rid edu];
rid_edu = rid_edu(rid_edu(:, 2)>0, :); % exclude -1 & -4
rid_edu = sortrows(rid_edu, 1);

% Append dummy columns
edu = cell(size(rid_edu, 1), 4);
edu(:, 1) = num2cell(rid_edu(:, 1));
edu(:, 2) = repmat({'bl'}, size(edu, 1), 1);
edu(:, 3) = repmat({'bl'}, size(edu, 1), 1);
edu(:, 4) = num2cell(rid_edu(:, 2));
