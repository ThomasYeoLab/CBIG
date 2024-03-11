function CBIG_pFIC_PNC_cognition_effect_prepare_input_wrapper(TC_path, subject_all, age_all, behavior_score_all)

% CBIG_pFIC_PNC_cognition_effect_prepare_input_wrapper(subject_all, age_all, behavior_score_all)
% This is a wrapper function to generate group-level FC and FCD for
% cognition effect analysis. Note that this function requires access to the subject's time
% series data, which are not open to public under the DUA. Thus this
% function is for CBIG internal usage only.
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - subject_all: the path to a file with 885 rows, each row
%   represents a subject ID sorted in ascending order according to age
%   - age_all: the path to a file with 885 rows, each row
%   represents a subject'age sorted in ascending order. 
%   - behavior_score_all: the path to a file with 885 rows, each row
%   represents a subject'overall accuracy score.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Note that the subjects should be sorted by age in an ascending order,
% the n-th row of age_all and behavior_score_all should
% CORRESPOND to the n-th subject from subject_all.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% Example:
% CBIG_pFIC_PNC_cognition_effect_prepare_input_wrapper(['/isilon/CSC1/Yeolab/Data/PNC/' ... 
%   'desikan_tc_from_vol/rest'], ...
%   '/home/shaoshi.z/storage/MFM/PNC/script/PNC_rest_subject_list_sorted_age.txt', ...
%   '/home/shaoshi.z/storage/MFM/PNC/script/rest_age_sort.txt', ...
%   '/home/shaoshi.z/storage/MFM/PNC/behav_score/rest/overall_acc.txt');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% read in the data
subject_all = dlmread(subject_all);
age_all = dlmread(age_all);
behavior_score_all = dlmread(behavior_score_all);

%% divide subjects into different age groups based on their age (12 months apart)
% make a temporary folder to store intermediate results
if exist('./tmp', 'dir')
    rmdir('tmp', 's');
    mkdir('tmp');
else
    mkdir('tmp');
end

counter = 1;
for i = min(age_all):12:max(age_all)
   	subject_sublist = [];
    for j = 1:length(subject_all)
        if age_all(j) >= i && age_all(j) < i + 12
            subject_sublist = [subject_sublist; subject_all(j)];
        end
    end
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_sublist, 'UniformOutput', false);
    fid = fopen(['tmp/subject_list_' num2str(counter) '.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    counter = counter + 1;
end

%% separate all the subjects into low and high-performance groups
% Within each age group (done by the pervious section), divide the subjects
% into a low and a high performance group based on the median of the
% cognitive score.

subject_info_low = [];
subject_info_high = [];

for i =  1:15
    subject_list = dlmread(['tmp/subject_list_' num2str(i) '.txt']); 
    subject_info_mat = zeros(length(subject_list), 3); % subject-id, behav, age
    for j = 1:length(subject_list)
        ind_high = find(subject_all == subject_list(j));
        subject_info_mat(j,1) = subject_list(j);
        subject_info_mat(j,2) = behavior_score_all(ind_high);
        subject_info_mat(j,3) = age_all(ind_high);
    end
    
    subject_info_mat_sorted = sortrows(subject_info_mat, 2);
    subject_low = subject_info_mat_sorted(1:ceil(size(subject_info_mat, 1)/2), :);
    subject_high = subject_info_mat_sorted(ceil(size(subject_info_mat, 1)/2)+1:end, :);
    
    subject_info_low = [subject_info_low; subject_low];
    subject_info_high = [subject_info_high; subject_high];
end
subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_high(:, 1), 'UniformOutput', false);
fid = fopen('tmp/subject_list_high.txt', 'w');
for j = 1:numel(subject_sublist_strings)
    fprintf(fid, '%s\n', subject_sublist_strings{j});
end
fclose(fid);

subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_low(:, 1), 'UniformOutput', false);
fid = fopen('tmp/subject_list_low.txt', 'w');
for j = 1:numel(subject_sublist_strings)
    fprintf(fid, '%s\n', subject_sublist_strings{j});
end
fclose(fid);

dlmwrite('tmp/overall_acc_high.txt', subject_info_high(:, 2), 'precision', 15);
dlmwrite('tmp/overall_acc_low.txt', subject_info_low(:, 2), 'precision', 15);

dlmwrite('tmp/age_high.txt', subject_info_high(:, 3), 'precision', 15);
dlmwrite('tmp/age_low.txt', subject_info_low(:, 3), 'precision', 15);

%% allocate space for each subgroup of low/high-performance groups
% High and low performance group each has 14 subgroups in TOTAL, the size of each
% subgroup varies a bit (either 31 or 32 subjects per subgroup)
low_performance_group_length = [32*ones(10, 1); 31*ones(4, 1)];
high_performance_group_length = [32*ones(7, 1); 31*ones(7, 1)];

%% generate subject lists for low/high-performance groups
ind_high = 1;
ind_low = 1;
for i = 1:14
    subject_info_sublist_high = subject_info_high(ind_high:ind_high+high_performance_group_length(i)-1, :);
    subject_info_sublist_low = subject_info_low(ind_low:ind_low+low_performance_group_length(i)-1, :);
    
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_sublist_high(:, 1), 'UniformOutput', false);
    fid = fopen(['tmp/subject_list_' num2str(i) '_high.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_sublist_low(:, 1), 'UniformOutput', false);
    fid = fopen(['tmp/subject_list_' num2str(i) '_low.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    
    dlmwrite(['tmp/age_' num2str(i) '_high.txt'], subject_info_sublist_high(:, 3), 'precision', 15);
    dlmwrite(['tmp/age_' num2str(i) '_low.txt'], subject_info_sublist_low(:, 3), 'precision', 15);
    dlmwrite(['tmp/overall_acc_' num2str(i) '_high.txt'], subject_info_sublist_high(:, 2), 'precision', 15);
    dlmwrite(['tmp/overall_acc_' num2str(i) '_low.txt'], subject_info_sublist_low(:, 2), 'precision', 15);
    ind_high = ind_high + high_performance_group_length(i);
    ind_low = ind_low + low_performance_group_length(i);
end

%% divide into training and validation set
if exist('../../replication/PNC/cognition_effect/input', 'dir')
   rmdir('../../replication/PNC/cognition_effect/input', 's'); 
   mkdir('../../replication/PNC/cognition_effect/input'); 
else
   mkdir('../../replication/PNC/cognition_effect/input'); 
end
mkdir('../../replication/PNC/cognition_effect/input/high_performance'); 
mkdir('../../replication/PNC/cognition_effect/input/low_performance'); 
for i = 1:14
    disp(['Preparing Input for Group ' num2str(i) '...'])
    mkdir(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) '/']); 
    mkdir(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) '/']); 
    
    subject_list_high = dlmread(['tmp/subject_list_' num2str(i) '_high.txt']);
    subject_list_low = dlmread(['tmp/subject_list_' num2str(i) '_low.txt']);
    age_high = dlmread(['tmp/age_' num2str(i) '_high.txt']);
    age_low = dlmread(['tmp/age_' num2str(i) '_low.txt']);
    behavior_score_high = dlmread(['tmp/overall_acc_' num2str(i) '_high.txt']);
    behavior_score_low = dlmread(['tmp/overall_acc_' num2str(i) '_low.txt']);
    
    rng(2000, 'twister')
    subject_info_high = [subject_list_high age_high behavior_score_high];
    subject_info_high_perm = subject_info_high(randperm(size(subject_info_high, 1)), :);
    subject_info_high_train = subject_info_high_perm(1:16, :);
    subject_info_high_validation = subject_info_high_perm(17:end, :);
    
    rng(2000, 'twister')
    subject_info_low = [subject_list_low age_low behavior_score_low];
    subject_info_low_perm = subject_info_low(randperm(size(subject_info_low, 1)), :);
    subject_info_low_train = subject_info_low_perm(1:16, :);
    subject_info_low_validation = subject_info_low_perm(17:end, :);
    
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_high_train(:, 1), 'UniformOutput', false);
    fid = fopen(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) ...
        '/training_subject_list.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_high_validation(:, 1), 'UniformOutput', false);
    fid = fopen(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) ...
        '/validation_subject_list.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    dlmwrite(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) ...
        '/training_subject_age.txt'], subject_info_high_train(:, 2), 'precision', 15);
    dlmwrite(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) ...
        '/validation_subject_age.txt'], subject_info_high_validation(:, 2), 'precision', 15);
    dlmwrite(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) ...
        '/training_subject_behav.txt'], subject_info_high_train(:, 3), 'precision', 15);
    dlmwrite(['../../replication/PNC/cognition_effect/input/high_performance/' num2str(i) ...
        '/validation_subject_bheav.txt'], subject_info_high_validation(:, 3), 'precision', 15);
    
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_low_train(:, 1), 'UniformOutput', false);
    fid = fopen(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) ...
        '/training_subject_list.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    subject_sublist_strings = arrayfun(@(x) num2str(x), subject_info_low_validation(:, 1), 'UniformOutput', false);
    fid = fopen(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) ...
        '/validation_subject_list.txt'], 'w');
    for j = 1:numel(subject_sublist_strings)
        fprintf(fid, '%s\n', subject_sublist_strings{j});
    end
    fclose(fid);
    dlmwrite(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) ...
        '/training_subject_age.txt'], subject_info_low_train(:, 2), 'precision', 15);
    dlmwrite(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) ...
        '/validation_subject_age.txt'], subject_info_low_validation(:, 2), 'precision', 15);
    dlmwrite(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) ...
        '/training_subject_behav.txt'], subject_info_low_train(:, 3), 'precision', 15);
    dlmwrite(['../../replication/PNC/cognition_effect/input/low_performance/' num2str(i) ...
        '/validation_subject_bheav.txt'], subject_info_low_validation(:, 3), 'precision', 15);
    
    % generate and save FC and FCD
    [training_FC, FCD] = CBIG_pFIC_generate_PNC_group_level_FC_FCD(TC_path, ['../../replication/PNC/' ... 
        'cognition_effect/input/high_performance/' num2str(i) '/training_subject_list.txt']);
    dlmwrite(['../../replication/PNC/cognition_effect/input/' ... 
        'high_performance/' num2str(i) '/FC_train.csv'], training_FC, 'delimiter', ',', 'precision', 15);
    save(['../../replication/PNC/cognition_effect/input/' ... 
        'high_performance/' num2str(i) '/FCD_train.mat'], 'FCD');
    [validation_FC, FCD] = CBIG_pFIC_generate_PNC_group_level_FC_FCD(TC_path, ['../../replication/PNC/' ... 
        'cognition_effect/input/high_performance/' num2str(i) '/validation_subject_list.txt']);
    dlmwrite(['../../replication/PNC/cognition_effect/input/' ... 
        'high_performance/' num2str(i) '/FC_validation.csv'], validation_FC, 'delimiter', ',', 'precision', 15);
    save(['../../replication/PNC/cognition_effect/input/' ... 
        'high_performance/' num2str(i) '/FCD_validation.mat'], 'FCD');
    
    [training_FC, FCD] = CBIG_pFIC_generate_PNC_group_level_FC_FCD(TC_path, ['../../replication/PNC/' ... 
        'cognition_effect/input/low_performance/' num2str(i) '/training_subject_list.txt']);
    dlmwrite(['../../replication/PNC/cognition_effect/input/' ... 
        'low_performance/' num2str(i) '/FC_train.csv'], training_FC, 'delimiter', ',', 'precision', 15);
    save(['../../replication/PNC/cognition_effect/input/' ... 
        'low_performance/' num2str(i) '/FCD_train.mat'], 'FCD');
    [validation_FC, FCD] = CBIG_pFIC_generate_PNC_group_level_FC_FCD(TC_path, ['../../replication/PNC/' ... 
        'cognition_effect/input/low_performance/' num2str(i) '/validation_subject_list.txt']);
    dlmwrite(['../../replication/PNC/cognition_effect/input/' ... 
        'low_performance/' num2str(i) '/FC_validation.csv'], validation_FC, 'delimiter', ',', 'precision', 15);
    save(['../../replication/PNC/cognition_effect/input/' ... 
        'low_performance/' num2str(i) '/FCD_validation.mat'], 'FCD');
end

rmdir('tmp', 's');
end