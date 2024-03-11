function CBIG_pFIC_GUSTO_prepare_input_wrapper(TC_path, sorted_subject_list_path, sorted_age_list_path, output_path)

% CBIG_pFIC_GUSTO_prepare_input_wrapper(TC_path, sorted_subject_list_path, sorted_age_list_path, output_path)
% This is a wrapper function to generate group-level FC and FCD for GUSTO
% cognitive effect analysis. Note that this function requires access to the subject's time
% series data, which are not open to public under the DUA. Thus this
% function is for CBIG internal usage only.
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - sorted_subject_list_path: the path to a file with 154 rows, each row
%   represents a subject ID sorted in ascending order according to age
%   - sorted_age_list_path: the path to a file with 154 rows, each row
%   represents a subject'age sorted in ascending order. This file has a row-to-row corresponance
%   with the above sorted subject list.
%   - output_path: the path to a folder where the generated subject lists
%   will be saved. If the folder does not exist, this script will attempt
%   to create that folder.
%
% Example:
% CBIG_pFIC_GUSTO_prepare_input_wrapper(['/isilon/CSC1/Yeolab/Data/GUSTO/pFIC/' ... 
%   'desikan_tc_noGSR/7.5'], ...
%   '/home/shaoshi.z/storage/MFM/GUSTO/behavior/subject_list/subject_list_sorted_age.txt', ...
%   '/home/shaoshi.z/storage/MFM/GUSTO/behavior/subject_list/age_list_sorted_age.txt', ...
%   '../../replication/GUSTO/input/'')
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% generate composite behavior scores using PCA
score_pca = CBIG_pFIC_PCA_GUSTO_behavior_score(sorted_subject_list_path);
[score_sort, ind] = sort(score_pca, 'descend');

%% divide subject list into low and high performance groups
fid = fopen(sorted_subject_list_path);
subject_list = textscan(fid, '%s');
subject_list = subject_list{1};
fclose(fid);
subject_list_sort = subject_list(ind);
age_list = dlmread(sorted_age_list_path);
age_list_sort = age_list(ind);

low_performance_subject_list = subject_list_sort(1:floor(length(subject_list)/2));
high_performance_subject_list = subject_list_sort(floor(length(subject_list)/2)+1:end);
low_performance_age_list = age_list_sort(1:floor(length(subject_list)/2));
high_performance_age_list = age_list_sort(floor(length(subject_list)/2)+1:end);
low_performance_behavior_list = score_sort(1:floor(length(subject_list)/2));
high_performance_behavior_list = score_sort(floor(length(subject_list)/2)+1:end);

%% generate training and validation split within low and high performance groups
rng(134, 'twister');

low_rand_ind = randperm(size(low_performance_subject_list, 1));
high_rand_ind = randperm(size(high_performance_subject_list, 1));

low_performance_subject_list_perm = low_performance_subject_list(low_rand_ind);
high_performance_subject_list_perm = high_performance_subject_list(high_rand_ind);
low_performance_age_list_perm = low_performance_age_list(low_rand_ind);
high_performance_age_list_perm = high_performance_age_list(high_rand_ind);
low_performance_behavior_list_perm = low_performance_behavior_list(low_rand_ind);
high_performance_behavior_list_perm = high_performance_behavior_list(high_rand_ind);

validation_low_subject_list = ...
    low_performance_subject_list_perm(ceil(length(low_performance_subject_list_perm)/2)+1:end);
validation_high_subject_list = ...
    high_performance_subject_list_perm(ceil(length(high_performance_subject_list_perm)/2)+1:end);
validation_low_age_list = ...
    low_performance_age_list_perm(ceil(length(low_performance_age_list_perm)/2)+1:end);
validation_high_age_list = ...
    high_performance_age_list_perm(ceil(length(high_performance_age_list_perm)/2)+1:end);
validation_low_behav_list = ...
    low_performance_behavior_list_perm(ceil(length(low_performance_behavior_list_perm)/2)+1:end);
validation_high_behav_list = ...
    high_performance_behavior_list_perm(ceil(length(high_performance_behavior_list_perm)/2)+1:end);

training_low_subject_list = low_performance_subject_list_perm(1:ceil(length(low_performance_subject_list_perm)/2));
training_high_subject_list = high_performance_subject_list_perm(1:ceil(length(high_performance_subject_list_perm)/2));
training_low_age_list = low_performance_age_list_perm(1:ceil(length(low_performance_age_list_perm)/2));
training_high_age_list = high_performance_age_list_perm(1:ceil(length(high_performance_subject_list_perm)/2));
training_low_behav_list = low_performance_behavior_list_perm(1:ceil(length(low_performance_behavior_list_perm)/2));
training_high_behav_list = high_performance_behavior_list_perm(1:ceil(length(high_performance_behavior_list_perm)/2));

%% prepare output folders
if ~exist([output_path '/high/'], 'dir')
    mkdir([output_path '/high/']);
end
if ~exist([output_path '/low/'], 'dir')
    mkdir([output_path '/low/']);
end

%% save subject list and corresponding age and behavior score
filePh = fopen([output_path '/high/training_subject_list.txt'],'w');
fprintf(filePh,'%s\n',training_high_subject_list{:});
fclose(filePh);
filePh = fopen([output_path '/high/validation_subject_list.txt'],'w');
fprintf(filePh,'%s\n',validation_high_subject_list{:});
fclose(filePh);
dlmwrite([output_path '/high/training_age_list.txt'], training_high_age_list);
dlmwrite([output_path '/high/validation_age_list.txt'], validation_high_age_list);
dlmwrite([output_path '/high/training_behav_list.txt'], training_high_behav_list);
dlmwrite([output_path '/high/validation_behav_list.txt'], validation_high_behav_list);

filePh = fopen([output_path '/low/training_subject_list.txt'],'w');
fprintf(filePh,'%s\n',training_low_subject_list{:});
fclose(filePh);
filePh = fopen([output_path '/low/validation_subject_list.txt'],'w');
fprintf(filePh,'%s\n',validation_low_subject_list{:});
fclose(filePh);
dlmwrite([output_path '/low/training_age_list.txt'], training_low_age_list);
dlmwrite([output_path '/low/validation_age_list.txt'], validation_low_age_list);
dlmwrite([output_path '/low/training_behav_list.txt'], training_low_behav_list);
dlmwrite([output_path '/low/validation_behav_list.txt'], validation_low_behav_list);

%% generate training and validation FC and FCD
[training_low_FC, FCD_train] = ...
    CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, [output_path '/low/training_subject_list.txt']);
[validation_low_FC, FCD_validation] = ...
    CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, [output_path '/low/validation_subject_list.txt']);
dlmwrite([output_path, '/low/FC_train.csv'], training_low_FC, 'delimiter', ',', 'precision', 15);
dlmwrite([output_path, '/low/FC_validation.csv'], validation_low_FC, 'delimiter', ',', 'precision', 15);
save([output_path, '/low/FCD_train.mat'], 'FCD_train');
save([output_path, '/low/FCD_validation.mat'], 'FCD_validation');

[training_high_FC, FCD_train] = ...
    CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, [output_path '/high/training_subject_list.txt']);
[validation_high_FC, FCD_validation] = ...
    CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, [output_path '/high/validation_subject_list.txt']);
dlmwrite([output_path, '/high/FC_train.csv'], training_high_FC, 'delimiter', ',', 'precision', 15);
dlmwrite([output_path, '/high/FC_validation.csv'], validation_high_FC, 'delimiter', ',', 'precision', 15);
save([output_path, '/high/FCD_train.mat'], 'FCD_train');
save([output_path, '/high/FCD_validation.mat'], 'FCD_validation');

end