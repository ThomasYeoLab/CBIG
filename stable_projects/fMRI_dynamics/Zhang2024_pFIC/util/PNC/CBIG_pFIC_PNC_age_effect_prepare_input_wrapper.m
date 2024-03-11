function CBIG_pFIC_PNC_age_effect_prepare_input_wrapper(sorted_subject_list_path, sorted_age_list_path, output_path)

% CBIG_pFIC_PNC_age_effect_prepare_input_wrapper(sorted_subject_list_path, sorted_age_list_path, output_path)
% This is a wrapper function to generate group-level FC and FCD for age
% effect analysis. Note that this function requires access to the subject's time
% series data, which are not open to public under the DUA. Thus this
% function is for CBIG internal usage only.
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
% CBIG_pFIC_PNC_age_effect_prepare_input_warapper( ...
% '/home/shaoshi.z/storage/MFM/PNC/script/PNC_rest_subject_list_sorted_age.txt', ...
% '/home/shaoshi.z/storage/MFM/PNC/script/rest_age_sort.txt', '../../replication/PNC/age_effect/input/')
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% split the subjects according to their age
CBIG_pFIC_PNC_age_effect_generate_subject_list(sorted_subject_list_path, ...
    sorted_age_list_path, [output_path, '/subject_list'])

%% generate 29 age groups
for group_number = 1:29
    disp(['Preparing input for PNC age effect group [' num2str(group_number) '] ...'])
    output_path_group = [output_path '/' num2str(group_number)];

    if ~exist(output_path_group, 'dir')
       mkdir(output_path_group)
    end

    subject_list = dlmread([output_path, '/subject_list/subject_list_' num2str(group_number) '.txt']);
    age_list = dlmread([output_path, '/subject_list/age_list_' num2str(group_number) '.txt']);
    subject_age_list = [subject_list age_list];
    rng(296, 'twister');
    subject_age_list_perm = subject_age_list(randperm(size(subject_age_list, 1)), :);
    subject_list_perm = subject_age_list_perm(:, 1);
    age_list_perm = subject_age_list_perm(:, 2);

    % 15 or 16 subjects in training set, always 15 subjects in validation set
    training_subject_list = subject_list_perm(16:end);
    training_subject_age = age_list_perm(16:end);
    validation_subject_list = subject_list_perm(1:15);
    validation_subject_age = age_list_perm(1:15);

    % save subject list and age list 
    dlmwrite([output_path_group '/training_subject_list.txt'], training_subject_list, 'precision', 20);
    dlmwrite([output_path_group '/validation_subject_list.txt'], validation_subject_list, 'precision', 20);

    dlmwrite([output_path_group '/training_subject_age.txt'], training_subject_age, 'precision', 20);
    dlmwrite([output_path_group '/validation_subject_age.txt'], validation_subject_age, 'precision', 20);

    % generate and save FC and FCD
    [training_FC, FCD] = CBIG_pFIC_generate_PNC_group_level_FC_FCD([output_path_group '/training_subject_list.txt']);
    dlmwrite([output_path_group, '/FC_train.csv'], training_FC, 'delimiter', ',', 'precision', 15);
    save([output_path_group, '/FCD_train.mat'], 'FCD');

    [validation_FC, FCD] = CBIG_pFIC_generate_PNC_group_level_FC_FCD( ...
        [output_path_group '/validation_subject_list.txt']);
    dlmwrite([output_path_group, '/FC_validation.csv'], validation_FC, 'delimiter', ',', 'precision', 15);
    save([output_path_group, '/FCD_validation.mat'], 'FCD');
end

end
