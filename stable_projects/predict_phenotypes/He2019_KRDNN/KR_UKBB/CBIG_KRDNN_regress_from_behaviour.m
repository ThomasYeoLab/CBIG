function CBIG_KRDNN_regress_from_behaviour(setup_file)

% CBIG_KRDNN_regress_from_behaviour(setup_file)
% 
% This function regress age, gender and motion from the behavioral measures.
% For age and gender, we did not regress out any thing.
%
% Inputs:
%   - setup_file
%     Path of the setup file (.txt file).
%     setup_file has the following fields:
%     1. Line 1
%        `data`, the location of your output and intermediate values. You can
%        also change this to any place you want.
%     2. Line 2
%        Path of subject list of UK biobank dataset. It should be a txt file
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
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[data_dir, subject_list, behavior_csv, HM_file, ~, ~, measure_set] = CBIG_KRDNN_mics_load_setup(setup_file);

mkdir(fullfile(data_dir, 'y', 'dis_y_regress'));
% load train/valid/test split
load(fullfile(data_dir, 'ukbb_subject_split.mat'))

measures_set = {measure_set};
c = 1;
y_regress_cell(:,1) = read_sub_list(subject_list);
for j=1:length(measures_set)
    curr_measure=measures_set{j};

    % extract behaviour measurements from restriced and unrestriced sheet
    measure_file = fullfile(data_dir, 'measures_lists', [curr_measure, '.txt']);
    measure_cell = read_sub_list(measure_file);
    [SubjectID_cell, Nusiance_num, measures_num] = read_measures_and_regressor(subject_list,...
        measure_cell, behavior_csv, HM_file);
   
    y_regress_cell(:,1) = SubjectID_cell;

    flag_regress = [0, 0, 1, 1];
    for m = 1:size(measures_num,2)
        c = c + 1;

        % do regression
        % train
        [measures_num_regress(subject_split.index_tra==1,m), beta{m}] = cv_regression_train(...
            measures_num(subject_split.index_tra==1,m), Nusiance_num(subject_split.index_tra==1,:), flag_regress(m));
        % valid
        measures_num_regress(subject_split.index_val==1,m) = cv_regression_test(...
            measures_num(subject_split.index_val==1,m), ...
            Nusiance_num(subject_split.index_val==1,:), beta{m}, flag_regress(m));
        % test
        measures_num_regress(subject_split.index_tes==1,m) = cv_regression_test(...
            measures_num(subject_split.index_tes==1,m), ...
            Nusiance_num(subject_split.index_tes==1,:), beta{m}, flag_regress(m));

        y_regress_cell(:,c) = num2cell(measures_num_regress(:,m));    
    end
    if(j == 1)
        measure_header = measure_cell;
    else
        measure_header = [measure_header measure_cell];
    end
    save(fullfile(data_dir, 'y', 'dis_y_regress', [curr_measure, '_y_regress.mat']),'measures_num_regress','beta');
end
end


function [y_hat, beta] = cv_regression_train(y,X,flag)
% y: Nx1
% X: NxP

% check y in case there is nan
nan_index_y = isnan(y);
% check X in case there is nan, exclude this subject for all covariates
nan_index_X = isnan(sum(X,2));

nan_index = (nan_index_y | nan_index_X);

% keep original y so that y has the same length for all measures
y_hat = y;
y_hat(nan_index,:) = NaN;

y(nan_index) = [];
X(nan_index,:) = [];

% demean nusiance regressors
X = bsxfun(@minus, X, mean(X));

% Add bias term for X

if flag == 1
    X = [ones(size(X,1),1) X];
else
    X = ones(size(X,1),1);
end

% y = X*beta + e
% 
% estimate beta
%  beta = ( X' * X )^-1 * X' * y
%  Px1  =  PxN  NxP      PxN  Nx1
beta = (X'*X)\(X'*y);
y_hat(~isnan(y_hat)) = y - X*beta;
end

function y_hat = cv_regression_test(y,X,beta,flag)
% y: Nx1
% X: NxP

% check y in case there is nan
nan_index_y = isnan(y);
% check X in case there is nan, exclude this subject for all covariates
nan_index_X = isnan(sum(X,2));

nan_index = (nan_index_y | nan_index_X);

% keep original y so that y has the same length for all measures
y_hat = y;
y_hat(nan_index,:) = NaN;

y(nan_index) = [];
X(nan_index,:) = [];

% demean nusiance regressors
X = bsxfun(@minus, X, mean(X));

% Add bias term for X

if flag == 1
    X = [ones(size(X,1),1) X];
else
    X = ones(size(X,1),1);
end

% y = X*beta + e
% 
y_hat(~isnan(y_hat)) = y - X*beta;
end

function [SubjectID_cell, regress_num, measures_num]=...
    read_measures_and_regressor(subject_list, measure_cell, behavior_csv, HM_file)

% [SubjectID_cell, regress_num, measures_num]=...
%   read_measures_and_regressor(subject_list, measure_cell, behavior_csv, HM_file)
% 
% This function regress age, gender and motion from the behavioral measures.
% For age and gender, we did not regress out any thing.
%
% Inputs:
%   - subject_list
%        Path of subject list of UK biobank dataset. It should be a txt file
%        that contains #subject of line (8868 in our paper), while each line
%        is the subject id of 1 subject.
%   - measure_cell
%        Cell of all behavioral measures to read from behavior_csv
%   - behavior_csv
%        Path of CSV file extracted from UK biobank dataset. You can
%        remove the unused measures from the csv to increase the loading
%        speed.
%   - HM_file
%        Path of the mean rfMRI head motion (data-field 25741) for each
%        subject. It should be a txt file that contains #subject of line,
%        while each line corresponds to the mean rfMRI head motion of a
%        subject. The order should follow the subject order of "subject_list"
%
% Outputs:
%   - SubjectID_cell
%        Cell of all subjects ID.
%   - regress_num
%        matrix of regress term for all subjects
%   - measures_num
%        matrix of behavioral measures term for all subjects
%


if(size(measure_cell,1)~=1)
    error('measure_cell should be a 1xK vector');
end

%% read in subject list
subj_list = read_sub_list(subject_list);

%% generate regressors
% age and gender
subjectname_str_cell = {'eid'};
nusiance_num_cell = {'age-2.0', '31-0.0'};
[SubjectID_cell, Age_Gender_num, ~, ~] = CBIG_parse_delimited_txtfile(...
    behavior_csv, subjectname_str_cell, nusiance_num_cell, 'eid', subj_list, ',');
Age_num = Age_Gender_num(:,1);
Gender_num = Age_Gender_num(:,2);

% FD
fileID = fopen(HM_file);
temp = textscan(fileID, '%f');
FD_num = temp{1};

regress_num = [Gender_num, Age_num, FD_num];

[~, measures_num, ~, ~] = CBIG_parse_delimited_txtfile(behavior_csv, [], measure_cell, 'eid', subj_list, ',');
end

function subj_list = read_sub_list(subject_text_list)
% this function will output a 1xN cell where N is the number of
% subjects in the text_list, each subject will be represented by one
% line in the text file
% NOTE: multiple runs of the same subject will still stay on the same
% line
% Each cell will contain the location of the subject, for e.g.
% '<full_path>/subject1_run1_bold.nii.gz <full_path>/subject1_run2_bold.nii.gz'
    fid = fopen(subject_text_list, 'r');
    i = 0;
    while(1);
        tmp = fgetl(fid);
        if(tmp == -1)
            break
        else
            i = i + 1;
            subj_list{i} = tmp;
        end
    end
    fclose(fid);
end
