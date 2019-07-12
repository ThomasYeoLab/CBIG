function [data_dir, subject_list_file, behavior_csv, HM_file, corr_mat_file, measure_dir, measure_set] ...
	= CBIG_KRDNN_mics_load_setup(setup_file)
% [data_dir, subject_list_file, behavior_csv, FD_txt_file, corr_mat_file, measure_dir, measure_set] ...
%	= CBIG_KRDNN_mics_load_setup(setup_file)
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
% Outputs:
%   - data_dir
%     Location to put output and intermediate values, See the description of
%     setup_file line 1
%   - subject_list_file
%     Path of subject list of UK biobank dataset. See the description of
%     setup_file line 2
%   - behavior_csv
%     Path of CSV file extracted from UK biobank dataset. See the description
%     of setup_file line 3
%   - HM_file
%     Path of the mean rfMRI head motion (data-field 25741) for each subject.
%     See the description of setup_file line 4
%   - corr_mat_file
%     Path of the functional connectivity matrix. See the description of
%     setup_file line 5
%   - measure_dir
%     Path of measure list text file
%   - measure_set
%     Name of the measures lists text file in measure_dir folder. See the
%     description of setup_file line 6
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


fileID = fopen(setup_file);
param = textscan(fileID, '%s');
fclose(fileID);

data_dir = param{1}{1};
subject_list_file = param{1}{2};
behavior_csv = param{1}{3};
HM_file = param{1}{4};
corr_mat_file = param{1}{5};
measure_dir = param{1}{6};
measure_set = param{1}{7};