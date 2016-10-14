function amyloid = CBIG_get_amyloid(rid_of_interest)

% amyloid = CBIG_get_amyloid(rid_of_interest)
% Fetch amyloid from .csv files
%
% Input:
%    - rid_of_interest: the RID's of subjects of interest
%
% Output:
%    - amyloid: 1st col is RID, 2nd col is (dummy) VISCODE, 3rd is (dummy)
%    EXAMDATE, 4th is value

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/UPENNBIOMK.csv', 'r');
data = textscan(file_id, ...
    '%*q %q %*q %q %*q %*q %*q %q %*[^\n]',...
    'delimiter', {',', '\n'}, 'HeaderLines', 1);
fclose(file_id);
data = [data{1, 1} data{1, 2} data{1, 3}];

% Remove rows with empty amyloid
keep_ind = ~cellfun(@isempty, data(:, 3));
data = data(keep_ind, :);

% Pick out RID
rid_col = cellfun(@str2num, data(:, 1));
keep_ind = ismember(rid_col, rid_of_interest);
data_rid = data(keep_ind, :);

% Add dummy columns
amyloid = cell(size(data_rid, 1), 4);
amyloid(:, 1) = num2cell(cellfun(@str2num, data_rid(:, 1)));
amyloid(:, 2) = data_rid(:, 2);
amyloid(:, 3) = repmat({''}, size(amyloid, 1), 1);
amyloid(:, 4) = num2cell(cellfun(@str2num, data_rid(:, 3)));
