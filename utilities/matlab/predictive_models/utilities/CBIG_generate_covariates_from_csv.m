function [covariates] = CBIG_generate_covariates_from_csv( csv_files, subject_header, ...
    covariate_names, covariate_types, subject_list, FD_file, DVARS_file, outname, delimiter )

% [covariates] = CBIG_generate_covariates_from_csv( csv_files, subject_header, ...
%     covariate_names, covariate_types, subject_list, FD_file, DVARS_file, outname, delimiter )
% 
% This function reads the covariates of all subjects in "subject_list" from
% multiple CSV files delimited by "delimiter", and combine with mean FD and
% DVARS of each subject (if necessary).
% 
% Inputs:
%   - csv_files
%     A string or cell arrays for multiple CSV files containing covariates,
%     e.g. age, gender, race. 
%     If "csv_files" is a string, it is the path to a list which contains
%     the full path to all csv files. Each line in this list corresponds to
%     one csv file. 
%     If "csv_files" are cell arrays, each cell array is a string
%     corresponding to the full path of a csv file.
% 
%   - subject_header
%     A string. The string should match with the header of the subject ID
%     column in the CSV files.
% 
%   - covariate_names
%     A string or cell arrays for multiple covariate names. Except for 'FD'
%     and 'DVARS', the covariate names should match with the headers names
%     in the CSV files.
%     If "covariate_names" is a string, it is the full path of a list. Each
%     line in this list corresponds to one covariate name.
%     If "covariate_names" are cell arrays, each cell array contains the
%     name of one covariate (string type).
% 
%   - covariate_types
%     A cell array of strings. The i-th array defines whether the i-th
%     covariate is 'categorical' or 'continuous'. Binary covariates are
%     also considered 'categorical'.
%     (1) If i-th covariate is binary (e.g. sex), a #subject x 1 vector
%     containing the values 0 or 1 will be generated and concatenated with
%     "covariates" matrix.
%     (2) If i-th covariate is categorical (with #UniqueValues > 2, e.g.
%     race), a #subjects x (#UniqueValues - 1) matrix will be generated and
%     concatenated with "covariates" matrix.
%     (3) If i-th covariate is continuous, the corresponding column will be
%     read and directly concatenated with "covariates" matrix.
%     
%   - subject_list
%     The full path of the subject list. Each line in this list is one
%     subject ID.
%   
%   - FD_file
%     If there is a need to regress 'FD' (framewise displacement) from the
%     target behavioral (or demographic) measures, y, the user must include
%     the covariate 'FD' in the 'covariate_list'. In this case, "FD_file"
%     is the full path of the mean framewise displacement (FD) of all
%     subjects. The number of lines in "FD_file" should be the same as the
%     number of lines in "subject_list". 
%     If the user does not need to regress 'FD' from y, then the input
%     variable 'FD_file' is not required and the user can pass in 'NONE' to
%     the function.
% 
%   - DVARS_file
%     If there is a need to regress 'DVARS' from the behavioral
%     (demographic) measures, y, the user must include the covariate
%     'DVARS' (or 'DV') in the 'covariate_list'. In this case, "DVARS_file"
%     is the full path of the mean DVARS of all subjects. The number of
%     lines in "DVARS_file" should be the same as the number of lines in
%     "subject_list".
%     If the user does not need to regress 'DVARS' from y, then the input
%     variable 'DVARS_file' is not required and the user can pass in 'NONE'
%     to the function.
% 
%   - outname
%     The full path of the output filename. If you don't want to save out
%     the covariates, pass in 'NONE'.
%    
%   - delimiter
%     The delimiter to separate columns in the CSV files (optional, default is ','). 
%     It can be a single character if all CSV files use the same delimiter. 
%     Or it can be a cell array whose length is the same as the number of
%     csv files, so that each CSV file corresponds to one delimiter.
% 
% Outputs:
%   - covariates
%     An #subjects x #covariates matrix with all covariates. The column
%     ordering follows the ordering in "covariate_names"
% 
% Example:
% [covariates] = CBIG_generate_covariates_from_csv( {'<path_to_1st_csv>/RESTRICTED.csv', ...
%     '<path_to_2nd_csv>/unrestricted.csv'}, 'Subject', {'Age_in_Yrs', 'Gender', 'FD', 'DVARS'}, ...
%     {'continuous', 'categorical', 'continuous', 'continuous'}, '<path_to_subject_list>/subject_953.txt', ...
%     '<path_to_FDfile>/FD_regressor_953.txt', '<path_to_DVARSfile>/DV_regressor_953.txt', ...
%     '<path_to_output>/covariates.mat', ',' )
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

if(ischar(covariate_names))
    [covariate_names, num_cov] = CBIG_text2cell(covariate_names);
elseif(iscellstr(covariate_names))
    num_cov = length(covariate_names);
else
    error('"covariate_names" should be a string or a cell array of strings.')
end

if(~iscellstr(covariate_types))
    error('"covariate_types" should be a cell array of strings.')
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

% For each covariate name, find which file (csv? FD? DVARS?) it corresponds to,
% and read that column
covariates = [];
for i = 1:num_cov
    cov_exist = 0;
    
    if(strcmp(covariate_names{i}, 'FD'))
        % If current covariate is FD, check if FD_file was passed in, and
        % check the length of FD
        cov_exist = 1;
        if(strcmpi(FD_file, 'none'))
            error('"FD" is one of the covariates. "FD_file" is necessary but not passed in.')
        else
            FD = dlmread(FD_file);
            if(length(FD)~=num_sub)
                error('Length of FD should equal to #subjects in "subject_list."')
            end
            covariates = [covariates FD];
        end
        
    elseif(strcmp(covariate_names{i}, 'DVARS') || strcmp(covariate_names{i}, 'DV'))
        % If current covariate is DVARS, check if DVARS_file was passed in,
        % and check the length of DVARS
        cov_exist = 1;
        if(strcmpi(DVARS_file, 'none'))
            error('"DVARS" is one of the covariates. "DVARS_file" is necessary but not passed in.')
        else
            DV = dlmread(DVARS_file);
            if(length(DV)~=num_sub)
                error('Length of DVARS should equal to #subjects in "subject_list."')
            end
            covariates = [covariates DV];
        end
        
    else
        % For each CSV file, check if current covariate name is part of its
        % headers
        curr_cov = CBIG_read_1header_from_Ncsv(csv_files, subject_header, subjects, ...
            covariate_names{i}, covariate_types{i}, delimiter);
        covariates = [covariates curr_cov];
        cov_exist = 1;
    end
    
    % If none of the above cases correspond to current covariate name,
    % throw an exception.
    if(cov_exist==0)
        error('Unknown covariate name: %s.', covariate_names{i});
    end
end

if(~strcmpi(outname, 'none'))
    outdir = fileparts(outname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir)
    end
    save(outname, 'covariates')
end

end

