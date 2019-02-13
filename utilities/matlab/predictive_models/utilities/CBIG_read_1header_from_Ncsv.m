function [data] = CBIG_read_1header_from_Ncsv( csv_files, subject_header, subjects, data_header, data_type, delimiter )

% [data] = CBIG_read_1header_from_Ncsv( csv_files, subject_header, subjects, data_header, data_type, delimiter )
% 
% Given a header name and multiple CSV files, this function find which CSV
% file this header belongs to, and read the corresponding column.
% 
% Possible application of this function includes reading a covariate or a
% behavioral measure from multiple CSV files. For example, in the Human
% Connectome Project, some behavioral measures are in the restricted csv
% file, some are in the unrestricted csv file. This function can be used
% for users without any knowledge about which CSV file that a behavioral
% measure was stored in.
% 
% Inputs:
%   - csv_files
%     Cell arrays. Each cell array contains a string corresponding to the
%     full path of a csv file.
% 
%   - subject_header
%     A string, the header of the subject ID column in CSV files.
% 
%   - subjects
%     Cell arrays. Each array is a string containing the ID of one subject.
% 
%   - data_header
%     A string, the header of the data column in CSV files.
% 
%   - data_type
%     Whether the data are classified as 'categorical' or 'continuous'.
%     Binary data are considered as 'categorical' type.
%     (1) If the data are binary, a #subject x 1 vector with 0 or 1 will
%     be generated as "data" output.
%     (2) If the data are categorical (with #UniqueValues > 2, e.g. race),
%     a #subjects x (#UniqueValues - 1) matrix will be generated as "data"
%     output.
%     (3) If the data are continuous, the corresponding column in the csv
%     file will be read as "data" output.
% 
%   - delimiter (optional)
%     Cell arrays whose length is the same as the number of csv files, such
%     that each array of delimiter can be applied to each corresponding CSV
%     file corresponds to one delimiter.
%     Default is ',' for every CSV file.
% 
% Outputs:
%   - data
%     A matrix or a vector depending on data_type.
%     If data_type = 'categorical' (binary), "data" will be a vector.
%     If data_type = 'categorical' (#UniqueValues > 2), "data" will be a
%     matrix with (#UniqueValues - 1) columns.
%     If data_type = 'continuous', "data" will be a vector.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_csv = length(csv_files);
num_sub = length(subjects);

if(~exist('delimiter', 'var'))
    delimiter(1:num_csv) = {','};
end

% Read headers of each csv file
headers = cell(num_csv, 1);
for c = 1:num_csv
    fid = fopen(csv_files{c}, 'r');
    headerline = fgetl(fid);
    fclose(fid);
    headers(c) = textscan(headerline, '%s', 'Delimiter', delimiter{c});
end

data = [];
for c = 1:num_csv
    if(any(strcmp(headers{c}, data_header)==1))
        if(strcmpi(data_type, 'categorical'))
            [curr_data] = CBIG_parse_delimited_txtfile(csv_files{c}, {data_header}, [], subject_header, subjects, delimiter{c});
            uniq_data = unique(curr_data);
            num_unq = length(uniq_data);
            
            if(num_unq > 2)
                % #UniqueValues > 2: create a binary vector for each unique
                % value (except for the last one to avoid collinearity)
                for u = 1:num_unq-1
                    curr_vec = strcmp(curr_data, uniq_data{u});
                    data = [data curr_vec];
                end
            elseif(num_unq == 2)
                % #UniqueValues == 2: create a single binary vector
                data = zeros(num_sub, 1);
                data(strcmp(curr_data, uniq_data{1})) = 1;
            else
                error('%s is same for all subjects.', data_header)
            end
        elseif(strcmpi(data_type, 'continuous'))
            [~, data] = CBIG_parse_delimited_txtfile(csv_files{c}, [], {data_header}, subject_header, subjects, delimiter{c});
        else
            error('Unknown data_type: %s. Choose from ''categorical'' or ''continuous''', data_type)
        end
        break
    end
end

if(isempty(data))
    error('Couldn''t find %s in any csv file. Please check the spelling.', data_header)
end


end