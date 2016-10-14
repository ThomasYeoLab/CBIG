function mmse = CBIG_get_mmse(rid_of_interest, predementedOnly)

% mmse = CBIG_get_mmse(rid_of_interest, predementedOnly)
% Fetch MMSE scores from .csv files
%
% Input:
%    - rid: the RID's of subjects of interest
%
% Output:
%    - mmse: 1st col is RID, 2nd is VISCODE, 3rd is EXAMDATE, 4th is score

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Read MMSE.csv
file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/MMSE.csv', 'r');
data = textscan(file_id, ...
    '%*q %*q %q %*q %*q %q %*q %*q %q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %q %*[^\n]',...
    'delimiter', {',', '\n'}, 'HeaderLines', 1); % read in RID, PTGENDER, and PTDOBYY
fclose(file_id);
data = [data{1, 1} data{1, 2} data{1, 3} data{1, 4}];
% Certain RID only
rid_col = cellfun(@str2num, data(:, 1));
keep_ind = ismember(rid_col, rid_of_interest);
data_rid = data(keep_ind, :);
% Remove those without a score
keep_ind = ones(size(data_rid, 1), 1);
for idx = 1:numel(keep_ind)
    if strcmp(data_rid{idx, 4}, '') % empty score
        keep_ind(idx) = 0;
    end
end
data_rid = data_rid(logical(keep_ind), :);

%-------------- Some subjects' VISCODE are missing
%%% Sol 1 - reference back to REGISTRY.csv
disp('Missing ADNI GO/2 dates -> reference back to REGISTRY.csv');
file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/REGISTRY.csv', 'r');
data = textscan(file_id, ...
    '%*q %*q %q %*q %*q %q %*q %*q %*q %*q %*q %q %*[^\n]', ...
    'delimiter', {',', '\n'}, 'HeaderLines', 1);
fclose(file_id);
data = [data{1, 1} data{1, 2} data{1, 3}];
for idx = 1:size(data_rid, 1)
    if strcmp(data_rid{idx, 3}, '') % empty EXAMDATE
        rid = data_rid{idx, 1};
        viscode = data_rid{idx, 2};
        missing_idx = logical(ismember(data(:, 1), rid)&ismember(data(:, 2), viscode));
        % Fill it up!
        data_rid{idx, 3} = data{missing_idx, 3};
    end
end
% Also unavailable in REGISTRY.csv? Discard!
keep_ind = ones(size(data_rid, 1), 1);
for idx = 1:numel(keep_ind)
    if strcmp(data_rid{idx, 3}, '') % empty examdate
        keep_ind(idx) = 0;
    end
end
data_rid = data_rid(logical(keep_ind), :);
% %%% Sol 2 - add an approximate date estimated from VISCODE
% disp('Missing ADNI GO/2 dates -> estimate from VISCODE');
% for idx = 1:size(data_rid, 1)
%     if strcmp(data_rid{idx, 3}, '') % empty EXAMDATE
%         rid  = data_rid{idx, 1};
%         sc_idx = logical(ismember(data_rid(:, 1), rid)&ismember(data_rid(:, 2), 'sc'));
%         sc_date_str = data_rid{sc_idx, 3};
%         sc_date = datenum(sc_date_str, 'yyyy-mm-dd');
%         viscode = data_rid{idx, 2};
%         months = strsplit(viscode, 'm');
%         days = str2double(months{1, 2})*30;
%         examdate = sc_date+days;
%         examdate_str = datestr(examdate, 'yyyy-mm-dd');
%         data_rid{idx, 3} = examdate_str;
%     end
% end
% %%% Sol 3 - just remove them
% disp('Missing ADNI GO/2 dates -> just discard');
% keep_ind = ones(size(data_rid, 1), 1);
% for idx = 1:numel(keep_ind)
%     if strcmp(data_rid{idx, 3}, '') % empty EXAMDATE
%         keep_ind(idx) = 0;
%     end
% end
% data_rid = data_rid(logical(keep_ind), :);

% Sort by RID
data_rid_sort = data_rid;
rid_col = cellfun(@str2num, data_rid(:, 1));
viscode_col = data_rid(:, 2);
examdate_col = data_rid(:, 3);
score_col = cellfun(@str2num, data_rid(:, 4));
[rid_col_sort, sort_ind] = sort(rid_col);
score_col_sort = score_col(sort_ind);
data_rid_sort(:, 1) = num2cell(rid_col_sort);
data_rid_sort(:, 2) = viscode_col(sort_ind);
data_rid_sort(:, 3) = examdate_col(sort_ind);
data_rid_sort(:, 4) = num2cell(score_col_sort);

% Remove -1 scores
data_rid_sort = data_rid_sort(score_col_sort>=0, :);

mmse = data_rid_sort;

% Only want pre-dementia scores
if predementedOnly == 1
    disp('Use only the pre-dementia data');
    predementia_ind = ones(size(mmse, 1), 1);
    rid_date_age = get_convDate_convAge(rid_of_interest);
    % Loop thru all convert dates
    for idx = 1:size(rid_date_age, 1)
        rid = rid_date_age{idx, 1};
        convDate = rid_date_age{idx, 2};
        % RID matched & date > convDate, already demented, so predementia_ind set to 0
        predementia_ind(ismember(cell2mat(mmse(:, 1)), rid)&(datenum(mmse(:, 3), 'yyyy-mm-dd')>=datenum(convDate))) = 0;
    end
    mmse = mmse(logical(predementia_ind), :);
% Only want post-dementia scores
% Truncate at conversion
elseif predementedOnly == -1
    disp('Use only the post-dementia data');
    postdementia_ind = ones(size(mmse, 1), 1);
    rid_date_age = get_convDate_convAge(rid_of_interest);
    % RID not even matched, never convereted, so postdementia_ind set to 0
    postdementia_ind(~ismember(cell2mat(mmse(:, 1)), cell2mat(rid_date_age(:, 1)))) = 0;
    % Loop thru all convert dates
    for idx = 1:size(rid_date_age, 1)
        rid = rid_date_age{idx, 1};
        convDate = rid_date_age{idx, 2};
        % RID matched & date < convDate, not demented yet, so postdementia_ind set to 0
        postdementia_ind(ismember(cell2mat(mmse(:, 1)), rid)&(datenum(mmse(:, 3), 'yyyy-mm-dd')<datenum(convDate))) = 0;
    end
    mmse = mmse(logical(postdementia_ind), :);
else
    disp('Use all the data we have');
end

end