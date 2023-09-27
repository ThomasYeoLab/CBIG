function [subdata_str_cell, subdata_num, headers, alldata_str] = CBIG_parse_delimited_txtfile(filename, fieldname_str_cell, fieldname_num_cell,filter_fieldname, filter_val_str_cell, delimiter)

% [subdata_str_cell, subdata_num, headers, alldata_str] = CBIG_parse_delimited_txtfile(filename, fieldname_str_cell, fieldname_num_cell,filter_fieldname, filter_val_str_cell, delimiter)
% Extract data from txtfile to matlab cell array or numerical matrix. This function assumes that the txtfile has a single headerline.
% And, you can choose the colume and row that you want to extract.
% 
%   [subdata_str_cell, subdata_num, headers, alldata_str] = CBIG_parse_delimited_txtfile(filename, fieldname_str_cell, fieldname_num_cell, ...
%     filter_fieldname, filter_val_str_cell, delimiter)
%   Input:
%       filename: file name
%       fieldname_str_cell  : the fieldname (column) in the header that you want to extract as a cell array
%       fieldname_num_cell  : the fieldname (column) in the header that you want to extract as a numerical matrix
%       filter_fieldname    : the fieldname (column) that you want to use as a filter to filter the row. e.g. SUBJECT_ID
%       filter_val_str_cell : the fieldname (row) that you want to extract based on the filter_fieldname column
%       delimiter           : delimiter (default ',')
%   Output:
%       subdata_str_cell    : cell array that you want
%       subdata_num         : numerical matrix that you want
%       headers             : header info
%       alldata_str         : all the data after filtering
%   Example: 
%       function [subdata_str_cell, subdata_num, headers, alldata_str] = CBIG_parse_delimited_txtfile(filename, fieldname_str_cell, fieldname_num_cell, ...
%    filter_fieldname, filter_val_str_cell, delimiter = '\t')
% 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if (nargin < 6)
    delimiter = ',';
end

if (nargin < 5)
    filter_fieldname = []; %grab all line
end

fid = fopen(filename, 'r');
%%Read the first line of the file containing the header information
headerline = fgetl(fid);
%close file
fclose(fid);
%Convert the column headers into individual elements stored in a cell array.
headers = textscan(headerline,'%s','Delimiter',delimiter);
tmp_str = [];

for i = 1:length(headers{1})
   
   tmp_str = [tmp_str '%s'];
    
end
fid = fopen(filename, 'r');
%%Read the first line of the file containing the header information
headerline = fgetl(fid);
tmp_alldata_str = textscan(fid,tmp_str,'Delimiter',delimiter);
fclose(fid);

NumSbj = size(tmp_alldata_str{1},1);

%first let's get the str data
tmp_subdata_str_cell = cell(NumSbj,length(fieldname_str_cell));
filter_col = 0;

if (~isempty(filter_fieldname))
    for ii = 1:length(headers{1})
        if (strcmp(headers{1}{ii}, filter_fieldname))
            filter_col = ii;
        break
        end
    end
    if (filter_col == 0)
        error(['Could not find column - ' filter_fieldname  ' -to filter rows by...']);
    else
        display(['Will filter rows base on ' filter_fieldname ' - column: ' num2str(filter_col) ' ...']);
    end
end


for c = 1:length(fieldname_str_cell)
    % find corresponding column
    cor_col = 0;
    for ii = 1:length(headers{1})
        if (strcmp(headers{1}{ii}, fieldname_str_cell{c}))
            cor_col = ii;
            break
        end
    end

    if (size(tmp_alldata_str{cor_col}, 1) ~= NumSbj)
        error('The read columns length is wrong!');
    end
    if (cor_col == 0)
        
        error(['Could not find field ' fieldname_str_cell{c}]);
    end

    for ii = 1:NumSbj
        tmp_subdata_str_cell{ii, c} = tmp_alldata_str{cor_col}{ii};
    end
end


%now let's get the numeric data
tmp_subdata_num = zeros(NumSbj,length(fieldname_num_cell));

for c = 1:length(fieldname_num_cell)
    % find corresponding column
    cor_col = 0;
    for ii = 1:length(headers{1})
        if (strcmp(headers{1}{ii}, fieldname_num_cell{c}))
            cor_col = ii;
            break
        end
    end

    if (size(tmp_alldata_str{cor_col}, 1) ~= NumSbj)
        error('The read columns length is wrong!');
    end

    if (cor_col == 0)
        
        error(['Could not find field ' fieldname_str_cell{c}]);
    end
    
    for ii = 1:NumSbj
        tmp_subdata_num(ii, c) = str2double(tmp_alldata_str{cor_col}{ii});
    end
end

if (filter_col > 0)
    
    include_line_ind = [];
    %%% let's determine the number of rows that pass the filter
    for i = 1:length(filter_val_str_cell)
        
        found_sbj = 0;
        for j = 1:NumSbj
            
            if (strcmp(filter_val_str_cell{i}, tmp_alldata_str{filter_col}{j}))
               include_line_ind = [include_line_ind, j];
               found_sbj = 1;
               break;
            end           
            
        end
        
        if (found_sbj == 0)
            display(['Could not find line corresponding to ' filter_val_str_cell{i}]);
        end        
    end
    
    display(['Number of lines remaining after filter operation : ' num2str(length(include_line_ind))]);
    
    subdata_str_cell = cell(length(include_line_ind),length(fieldname_str_cell));
    subdata_num = zeros(length(include_line_ind),length(fieldname_num_cell));
    alldata_str = cell(size(tmp_alldata_str));
    
    for ii = 1:length(include_line_ind)
        for jj = 1:length(fieldname_str_cell)
       
            subdata_str_cell{ii,jj} = tmp_subdata_str_cell{include_line_ind(ii), jj};
                        
        end
    end
    
    for ii = 1:length(include_line_ind)
        for jj = 1:length(fieldname_num_cell)
       
            subdata_num(ii,jj) = tmp_subdata_num(include_line_ind(ii), jj);
                        
        end
    end
    for ii = 1:length(tmp_alldata_str)
       
        for jj = 1:length(include_line_ind)
            if (include_line_ind(jj) <= length(tmp_alldata_str{ii}))
                alldata_str{ii}{jj} = tmp_alldata_str{ii}{include_line_ind(jj)};        
            else
                display(['There may be a mismatch for the column corresponding to ' headers{1}{ii}]);
                display(['No worries if you have not explicitly asked for this field']);
                alldata_str{ii}{jj} = '';
            end
        end
    end
    
   
else
    subdata_str_cell = tmp_subdata_str_cell;
    subdata_num = tmp_subdata_num;
    alldata_str = tmp_alldata_str;
end

return
