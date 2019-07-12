function CBIG_KRDNN_generate_K_test(setup_file)

% CBIG_KRDNN_generate_K_test(setup_file)
% 
% This function generate kernel for testing set for kernel ridge regression
% algorithm.
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

[data_dir, ~, ~, ~, corr_mat_file, ~, ~] = CBIG_KRDNN_mics_load_setup(setup_file);

FSM_set = {'corr'};

%% load correlation matrix of full set
load(corr_mat_file)

% load train/valid/test split
load(fullfile(data_dir, 'ukbb_subject_split.mat'))

% only select the lower triangular
ltri_index=find(tril(ones(size(corr_mat,1),size(corr_mat,2)),-1)==1);
corr_mat_KbyS = reshape(corr_mat, size(corr_mat,1)*size(corr_mat,2), size(corr_mat,3));
corr_mat_KbyS = corr_mat_KbyS(ltri_index,:);
num_subject = size(corr_mat_KbyS,2);

corr_mat_KbyS_test = corr_mat_KbyS(:,subject_split.index_tes==1);
corr_mat_KbyS_train = corr_mat_KbyS(:,subject_split.index_tra==1);

for i=1:length(FSM_set)
    curr_FSM = FSM_set{1,i};
    FSM_out = fullfile(data_dir, 'FSM_cv', 'dis_test');
    mkdir(FSM_out);
    FSM_dir = fullfile(FSM_out, ['FSM_' FSM_set{1,i} '.mat']);
    if (~exist(FSM_dir,'file'))
        if(strcmp(curr_FSM,'corr'))
            curr_FSM_type = curr_FSM;
            scale = [];
        else
            index = regexp(curr_FSM,'_');
            scale = str2double(curr_FSM(index+1:end));
            curr_FSM_type = curr_FSM(1:index-1);
        end
        fprintf('current FSM type is %s \n',curr_FSM);
        FSM=CBIG_KRDNN_kernel_with_scale( ...
            corr_mat_KbyS_train,corr_mat_KbyS_test, ...
            subject_split.index_tra==1,subject_split.index_tes==1, ...
            curr_FSM_type,scale);

        save(FSM_dir,'FSM');
    end
end
