function CBIG_pFIC_PNC_age_effect_generate_subject_list(sorted_subject_list_path, sorted_age_list_path, output_path)

% CBIG_pFIC_PNC_age_effect_generate_subject_list(sorted_subject_list_path, sorted_age_list_path, output_path)
% This is generate 29 lists of subjects that will be used for PNC age
% effect analysis. 885 subjects are first sorted according to their age in
% an ascending order. Then the 885 subjects are split into 29 lists of 30
% or 31 subjects each. 
% Input:
%   - sorted_subject_list_path: the path to a file with 885 rows, each row
%   represents a subject ID sorted in ascending order according to age
%   - sorted_age_list_path: the path to a file with 885 rows, each row
%   represents a subject'age sorted in ascending order. This file has a row-to-row corresponance
%   with the above sorted subject list.
%   - output_path: the path to a folder where the generated subject lists
%   will be saved. If the folder does not exist, this script will attempt
%   to create that folder.
%
% Example:
% CBIG_pFIC_PNC_age_effect_generate_subject_list( ...
% '/home/shaoshi.z/storage/MFM/PNC/script/PNC_rest_subject_list_sorted_age.txt', ...
% '/home/shaoshi.z/storage/MFM/PNC/script/rest_age_sort.txt', '../../replication/PNC/age_effect/input/subject_list/')
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% read in all the subject IDs
fid = fopen(sorted_subject_list_path, 'r');
stringPattern = textscan(fid,'%s');
fclose(fid);
subject_list_all_subject = string(stringPattern{:});

%% read in all the subject age
fid = fopen(sorted_age_list_path, 'r');
stringPattern = textscan(fid,'%s');
fclose(fid);
age_list_all_subject = string(stringPattern{:});

%% generate 29 groups of subjects, 31 subjects in the first 15 groups, 30 in the remaining 14 groups
subject_list_length = [31*ones(15, 1); 30*ones(14, 1)];

if ~exist(output_path, 'dir')
   mkdir(output_path) 
end

counter = 1;
for i = 1:length(subject_list_length)
    subject_list = subject_list_all_subject(counter:counter+subject_list_length(i)-1);
    age_list = age_list_all_subject(counter:counter+subject_list_length(i)-1);
    fid = fopen([output_path, '/subject_list_' num2str(i) '.txt'], 'w');
    fprintf(fid, '%s\n',subject_list{:});
    fclose(fid);
    fid = fopen([output_path, '/age_list_' num2str(i) '.txt'], 'w');
    fprintf(fid, '%s\n',age_list{:});
    fclose(fid);
    counter = counter + subject_list_length(i);
end

end
