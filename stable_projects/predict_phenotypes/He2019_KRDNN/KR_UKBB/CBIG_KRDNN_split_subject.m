function CBIG_KRDNN_split_subject(setup_file, rng_num, num_val_tes)

% CBIG_KRDNN_split_subject(setup_file, rng_num, num_val_tes)
% 
% This function splits subjects into training, validation and testing set.
%
% Inputs:
%   - setup_file
%     Path of the setup file (.txt file).
%     setup_file has the following fields:
%     1. Line 1
%        `data`, the location of your output and intermediate values. You can
%        also change this to any place you want.
%     2. Line 2
%        Path of subject list of setup_file dataset. It should be a txt file
%        that contains #subject of line (8868 in our paper), while each line
%        is the subject id of 1 subject.
%     3. Line 3
%        Path of CSV file extracted from UK biobank dataset. You can
%        remove the unused measures from the csv to increase the loading
%        speed.
%     4. Line 4
%        Path of the mean rfMRI head motion (data-field 25741) for each
%        subject. It should be a txt file that contains #subject of line,
%        while each line corresponds to the mean rfMRI head motion of a
%        subject. The order should follow the subject order of "subject_list"
%     5. Line 5
%        Path of the functional connectivity matrix. A matrix "corr_mat"
%        is assumed to be saved in this file. "corr_mat" should be a 3-D
%        matrix with dimension of #ROI x #ROI x #subjects. Since it is a
%        connectivity matrix, only the lower-triangular off-diagonal entries
%        will be considered as features because the connectivity matrix is
%        symmetric. The order of "corr_mat" should follow the subject order
%        of "subject_list" 
%     6. Line 6
%        Path of measure list text file
%     7. Line 7
%        Name of the measures lists text file in measure_dir folder. You
%        do not need to change this unless you want to change the measure list
%   - rng_num
%     Number (integer) of random number generator
%   - num_val_tes
%     Number (integer) of subjects for validation and testing set.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% process num_val_tes
if ~isnumeric(num_val_tes)
    num_val_tes = str2double(num_val_tes);
end
if ~isscalar(num_val_tes)
    num_val_tes = num_val_tes(1);
end
if floor(num_val_tes) ~= num_val_tes
    num_val_tes = floor(num_val_tes);
end

% process rng_num
if ~isnumeric(rng_num)
    rng_num = str2double(rng_num);
end
if ~isscalar(rng_num)
    rng_num = rng_num(1);
end
if floor(rng_num) ~= rng_num
    rng_num = floor(rng_num);
end
rng(rng_num, 'twister');

[data_dir, subject_list_file, behavior_csv, ~, ~, measure_dir, ~] = CBIG_KRDNN_mics_load_setup(setup_file);

% create data_dir
mkdir(data_dir)

% copy measure list to data_dir
measure_list_orig = fullfile(measure_dir, 'measures_lists');
measure_list_new = fullfile(data_dir, 'measures_lists');
copyfile(measure_list_orig, measure_list_new)

% read out subject file
fileID = fopen(subject_list_file);
temp = textscan(fileID, '%d');
subject_list = temp{1};

% Train/Validation/Test split
num_total = size(subject_list, 1);
num_valid = num_val_tes;
num_test = num_val_tes;
num_train = num_total - num_test - num_valid;

% generate list
subject_list_rand = subject_list(randperm(num_total));
subject_list_valid = sort(subject_list_rand(1:num_valid));
subject_list_test = sort(subject_list_rand(num_valid+1:num_valid+num_test));
subject_list_train = sort(subject_list_rand(end-num_train+1:end));

% generate index list
subject_index_train = zeros(num_total,1);
for i = 1:size(subject_list_train,1)
    subject_index_train = subject_index_train + (subject_list_train(i) == subject_list);
end
subject_index_train = logical(subject_index_train);

subject_index_valid = zeros(num_total,1);
for i = 1:size(subject_list_valid,1)
    subject_index_valid = subject_index_valid + (subject_list_valid(i) == subject_list);
end
subject_index_valid = logical(subject_index_valid);

subject_index_test = zeros(num_total,1);
for i = 1:size(subject_list_test,1)
    subject_index_test = subject_index_test + (subject_list_test(i) == subject_list);
end
subject_index_test = logical(subject_index_test);

% save them in strucuture
subject_split = struct('list_tra', subject_list_train, 'list_val', subject_list_valid,...
    'list_tes', subject_list_test, 'index_tra', subject_index_train,...
    'index_val', subject_index_valid, 'index_tes', subject_index_test);

% save out
save(fullfile(data_dir, 'ukbb_subject_split.mat'), 'subject_split');

end
