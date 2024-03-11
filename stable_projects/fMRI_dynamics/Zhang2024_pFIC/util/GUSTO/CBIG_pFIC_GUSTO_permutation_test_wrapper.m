function permutation_test_input_struct = CBIG_pFIC_GUSTO_permutation_test_wrapper(TC_path, sorted_subject_list_path)

% permutation_test_input_struct = CBIG_pFIC_GUSTO_permutation_test_wrapper(TC_path, sorted_subject_list_path)
% This function generates 100 sets of FCs and FCD CDFs for permutation test. Note that this function requires 
% access to the subject's time series data, which are not open to public under the DUA. Thus this
% function is for CBIG internal usage only.
%
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - sorted_subject_list_path: the path to a text file (.txt) with 154 rows, each row
%   represents a subject ID sorted in ascending order according to age
% Output:
%   - permutation_test_input_struct: a structure with 8 fields, 4 fields
%   of FCs and 4 fields of FCD CDFs. Both FCs and FCD CDFs are further
%   divided into training and validation set of low and high performance groups.
%   The index of the last dimension of FC and FCD CDF (3rd
%   dimension for FC, 2nd dimension for FCD) establish the correspondance between the 2.
%
% Example:
% permutation_test_input_struct = CBIG_pFIC_GUSTO_permutation_test_wrapper(['/isilon/CSC1/Yeolab/Data/GUSTO/pFIC/' ... 
%   'desikan_tc_noGSR/7.5'], '/home/shaoshi.z/storage/MFM/GUSTO/behavior/subject_list/subject_list_sorted_age.txt');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% initialization
FC_train_high = zeros(68, 68, 100);
FC_validation_high = zeros(68, 68, 100);
FC_train_low = zeros(68, 68, 100);
FC_validation_low = zeros(68, 68, 100);

FCD_CDF_train_high = zeros(10000, 100);
FCD_CDF_validation_high = zeros(10000, 100);
FCD_CDF_train_low = zeros(10000, 100);
FCD_CDF_validation_low = zeros(10000, 100);

%% load subject list
fid = fopen(sorted_subject_list_path);
subject_list = textscan(fid, "%s");
subject_list = subject_list{1};
fclose(fid);

if exist('./tmp', 'dir')
    rmdir('tmp', 's');
    mkdir('tmp');
else
    mkdir('tmp');
end

%% generate 100 sets permutation test input FCs and FCDs
for seed = 1:100
    disp(['Preparing input FC and FCD CDF for permutation test [' num2str(seed) '] ...'])
    rng(seed, 'twister');
    ind = randperm(length(subject_list));
    subject_list_perm = subject_list(ind);
    
    low_subject_list = subject_list_perm(1:floor(length(subject_list)/2));
    high_subject_list = subject_list_perm(floor(length(subject_list)/2)+1:end);  
%% generate training and validation FC and FCD
    rng(1, 'twister');
    low_rand_ind = randperm(size(low_subject_list, 1));
    high_rand_ind = randperm(size(high_subject_list, 1));

    low_subject_list_perm = low_subject_list(low_rand_ind);
    high_subject_list_perm = high_subject_list(high_rand_ind);
    validation_low_subject_list = low_subject_list_perm(ceil(length(low_subject_list_perm)/2)+1:end);
    validation_high_subject_list = high_subject_list_perm(ceil(length(high_subject_list_perm)/2)+1:end);

    training_low_subject_list = low_subject_list_perm(1:ceil(length(low_subject_list_perm)/2));
    training_high_subject_list = high_subject_list_perm(1:ceil(length(high_subject_list_perm)/2));
    filePh = fopen('tmp/high_training_subject_list.txt','w');
    fprintf(filePh,'%s\n',training_high_subject_list{:});
    fclose(filePh);
    filePh = fopen('tmp/high_validation_subject_list.txt','w');
    fprintf(filePh,'%s\n',validation_high_subject_list{:});
    fclose(filePh);

    filePh = fopen('tmp/low_training_subject_list.txt','w');
    fprintf(filePh,'%s\n',training_low_subject_list{:});
    fclose(filePh);
    filePh = fopen('tmp/low_validation_subject_list.txt','w');
    fprintf(filePh,'%s\n',validation_low_subject_list{:});
    fclose(filePh);
    [FC_train_high(:, :, seed), FCD_CDF_train_high(:, seed)] = ...
        CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, 'tmp/high_training_subject_list.txt');
    [FC_train_low(:, :, seed), FCD_CDF_train_low(:, seed)] = ...
        CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, 'tmp/low_training_subject_list.txt');
    [FC_validation_high(:, :, seed), FCD_CDF_validation_high(:, seed)] = ...
            CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, 'tmp/high_validation_subject_list.txt');
    [FC_validation_low(:, :, seed), FCD_CDF_validation_low(:, seed)] = ...
            CBIG_pFIC_generate_GUSTO_group_level_FC_FCD(TC_path, 'tmp/low_validation_subject_list.txt');
end

%% prepare output
permutation_test_input_struct.FC_train_high = FC_train_high;
permutation_test_input_struct.FC_validation_high = FC_validation_high;
permutation_test_input_struct.FC_train_low = FC_train_low;
permutation_test_input_struct.FC_validation_low = FC_validation_low;

permutation_test_input_struct.FCD_CDF_train_high = FCD_CDF_train_high;
permutation_test_input_struct.FCD_CDF_validation_high = FCD_CDF_validation_high;
permutation_test_input_struct.FCD_CDF_train_low = FCD_CDF_train_low;
permutation_test_input_struct.FCD_CDF_validation_low = FCD_CDF_validation_low;

rmdir('tmp', 's');
end


