function [rid_conv, rid_nonconv] = CBIG_find_conv_nonconv(rid_of_interest, years)

% [rid_conv, rid_nonconv] = CBIG_find_conv_nonconv(rid_of_interest, years)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

file_id = fopen('//share/users/imganalysis/yeolab/data/ADNI_mert/scripts/DXSUM_PDXCONV_ADNIALL.csv', 'r');


data = textscan(file_id, ...
    '%*q %*q %q %*q %*q %q %*q %*q %q %q %q %*[^\n]',...
    'delimiter', {',', '\n'}, 'HeaderLines', 1); % read in RID, PTGENDER, and PTDOBYY
fclose(file_id);
data = [data{1, 1} data{1, 2} data{1, 3} data{1, 4} data{1, 5}];

% Certain RID only
rid_col = cellfun(@str2num, data(:, 1));
keep_ind = ismember(rid_col, rid_of_interest);
data_part = data(keep_ind, :);

% Sort by RID
rid_col = cellfun(@str2num, data_part(:, 1));
[~, sort_ind] = sort(rid_col); % sort
data_part_sort = data_part(sort_ind, :);

% %-------------- Convert in all data; nonconvert in all data & data > years
% 
% % Loop through subjects to see convert or not
% rid_conv = [];
% rid_nonconv = [];
% for idx = 1:numel(rid_of_interest)
%     rid = rid_of_interest(idx);
%     keep_ind = ismember(data_part_sort(:, 1), num2str(rid));
%     dates = data_part_sort(keep_ind, 3);
%     dxchange = str2num(cell2mat(data_part_sort(keep_ind, 4))); % ADNI-1 convention
%     dxcurren = str2num(cell2mat(data_part_sort(keep_ind, 5))); % ADNI-GO/2 convention
%     % In ANDI-1 window, already becomes AD; or
%     % In ADNI-GO/2 window, 3=Stable:AD-AD or 5=Conv:MCI-AD
%     if sum(dxcurren==3)~=0 || sum(dxchange==3|dxchange==5)~=0 % if convert
%         rid_conv = [rid_conv; rid];
%     elseif (datenum(dates{end})-datenum(dates{1}))/365 > years % if not convert within years
%         rid_nonconv = [rid_nonconv; rid];
%     end
% end

%-------------- Convert within years; nonconvert within years, data > years

% Loop through subjects to see convert or not
rid_conv = [];
rid_nonconv = [];
for idx = 1:numel(rid_of_interest)
    rid = rid_of_interest(idx);
    keep_ind = ismember(data_part_sort(:, 1), num2str(rid));
    dates = data_part_sort(keep_ind, 3);
    if (datenum(dates{end})-datenum(dates{1}))/365 > years % data > years
        keep_ind = find(keep_ind);
        % Remove data until years
        while (datenum(dates{end})-datenum(dates{1}))/365 > years+0.5
            dates(end) = [];
            keep_ind(end) = [];
        end
        tmp = zeros(size(data_part_sort, 1), 1);
        tmp(keep_ind) = 1;
        keep_ind = logical(tmp);
        dxchange = str2num(cell2mat(data_part_sort(keep_ind, 4))); % ADNI-1 convention
        dxcurren = str2num(cell2mat(data_part_sort(keep_ind, 5))); % ADNI-GO/2 convention
        % In ANDI-1 window, already becomes AD; or
        % In ADNI-GO/2 window, 3=Stable:AD-AD or 5=Conv:MCI-AD
        if sum(dxcurren==3)~=0 || sum(dxchange==3|dxchange==5)~=0 % if convert
            rid_conv = [rid_conv; rid];
        else % if not
            rid_nonconv = [rid_nonconv; rid];
        end
    else % data < years
        dxchange = str2num(cell2mat(data_part_sort(keep_ind, 4))); % ADNI-1 convention
        dxcurren = str2num(cell2mat(data_part_sort(keep_ind, 5))); % ADNI-GO/2 convention
        if sum(dxcurren==3)~=0 || sum(dxchange==3|dxchange==5)~=0 % if convert
            rid_conv = [rid_conv; rid];
        end
    end
    
end
