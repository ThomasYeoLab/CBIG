function [y] = CBIG_read_y_from_csv( csv_files, subject_header, y_names, y_types, subject_list, outname, delimiter )

% [y] = CBIG_read_y_from_csv( csv_files, subject_header, y_names, y_types, subject_list, outname, delimiter )
% 
% This function reads y (i.e. the target variables to predict, e.g.
% behaviors, disease groups) for any predictive model from multiple
% "csv_files" of all  subject in "subject_list". 
% 
% Inputs:
%   - csv_files
%     A string or cell arrays for multiple CSV files containing the target
%     variables to predict (e.g. behaviors, disease groups, ...). 
%     If the user has a text list which stores the path to all the csv
%     files, then "csv_files" can be the full path to this text list passed
%     in to the function as a string.
%     If "csv_files" are cell arrays, each cell array would contain a
%     string corresponding to the full path of a csv file.
% 
%   - subject_header
%     A string, the header of the subject ID column in CSV files.
% 
%   - y_names
%     A string or cell arrays containing the names of the multiple target
%     variables to be predicted (e.g. PMAT24_A_CR (fluid intelligence)).
%     The names should correspond to some header names in the CSV files.
%     If "y_names" is a string, then it should contain the full path to a
%     text list and each line in this list should correspond to a target
%     variable name.
%     If "y_names" are cell arrays, each cell array should contain one
%     variable name (string).
% 
%   - y_types
%     A  cell array of strings. The i-th array defines whether the i-th
%     target variable to predict is 'categorical' or 'continuous'. Binary
%     variables are considered as 'categorical'.
%     (1) If i-th target variable is binary (e.g. sex), a #subject x 1
%     vector with 0 or 1 will be generated and concatenated with the output
%     "y" matrix.
%     (2) If i-th target variable is categorical (with #UniqueValues > 2,
%     e.g. Alzhemier's diseas, MCI, control), a #subjects x 1 vector with
%     values from 1 to #UniqueValues will be generated and concatenated
%     with the output "y"  matrix.
%     (3) If i-th variable is continuous, the corresponding column will be
%     read and directly concatenated with "y" matrix.
%     
%   - subject_list
%     The full path of the subject list. Each line in this list corresponds
%     to one subject ID.
% 
%   - outname
%     The full path of the output filename. If you don't want to save out
%     the output matrix y, pass in 'NONE'.
%    
%   - delimiter
%     The delimiter to separate columns in the CSV files (optional, default is ','). 
%     It can be a single character if all CSV files use the same delimiter. 
%     Or it can be a cell array whose length is the same as the number of
%     csv files, such that each delimiter can be applied to its
%     corresponding CSV file.
% 
% Outputs:
%   - y
%     An #subjects x #TargetVariables matrix with all the target variables
%     to be predicted.
% 
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% check input format
if(ischar(csv_files))
    [csv_files, num_csv] = CBIG_text2cell(csv_files);
elseif(iscellstr(csv_files))
    num_csv = length(csv_files);
else
    error('"csv_files" should be a string or a cell array of strings.')
end

if(ischar(y_names))
    [y_names, num_y] = CBIG_text2cell(y_names);
elseif(iscellstr(y_names))
    num_y = length(y_names);
else
    error('"y_names" should be a string or a cell array of strings.')
end

if(~iscellstr(y_types))
    error('"y_types" should be a cell array of strings.')
end

if(~exist('delimiter', 'var'))
    delimiter(1:num_csv) = {','};
elseif(ischar(delimiter))
    delim_temp = delimiter;
    clear delimiter
    delimiter(1:num_csv) = {delim_temp};
elseif(~iscellstr(delimiter))
    error('"delimiter" should be a character or a cell array of strings.')
end

% Read subject IDs
[subjects, num_sub] = CBIG_text2cell(subject_list);

% Read y
y = [];
for i = 1:num_y
    curr_y = CBIG_read_1header_from_Ncsv( csv_files, subject_header, subjects, y_names{i}, y_types{i}, delimiter );
    
    % if y is categorical with more than 2 unique values, say P unique values
    % curr_y will be a #subjects x (P-1) matrix.
    % Hence all the columns need to be combined.
    if(size(curr_y,2) > 1)
        curr_y_temp = curr_y; 
        curr_y = zeros(num_sub, 1);
        % sum up (1st column x 1), (2nd column x 2), ..., ((P-1)-th column x (P-1))
        for col = 1:size(curr_y_temp, 2)
            curr_y = curr_y + curr_y_temp(:,col) .* col;
        end
        % and (indices not assigned to any of the P-1 columns x P)
        curr_y = curr_y + (1 - sum(curr_y_temp,2)) .* (size(curr_y_temp, 2) + 1);
    end
    y = [y curr_y];
end

if(~strcmpi(outname, 'none'))
    outdir = fileparts(outname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir)
    end
    save(outname, 'y')
end

end

