function CBIG_MM_HCP_save_data(out_dir, HCP_dir, measure_list_dir)

% CBIG_MM_HCP_save_data(List, out_dir, roi_400_nii)
% 
% This function save all data needed for HCP S1200
%
% Inputs:
%   - out_dir
%     Path of the your output data. You can also change this
%     to any place you want.
%
%   - HCP_dir
%     Path of the phenotypes csv location.
%
%   - measure_list_dir
%     Path of the HCP S1200 phenotypes list directory.
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

unrestricted_csv = fullfile(HCP_dir, 'subject_measures',...
    'unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv');
restricted_csv = fullfile(HCP_dir, ['restricted_hc' 'p_data'],...
    'RESTRICTED_jingweili_4_12_2017_1200subjects.csv');
subject_list = fullfile(out_dir, 'subject_list_FC_S1200_1094_210120.txt');

% get matrix of behavioral measures
csv_files = {unrestricted_csv};
subject_header = 'Subject';
measures_set = {'Cognitive','Personality_Task','Social_Emotion'};
y_names = {};
for i = measures_set
    name = i{1};
    measure_list = fullfile(measure_list_dir, [name '_unrestricted.txt']);
    temp = read_sub_list(measure_list);
    y_names = [y_names, temp];
end
y_names(ismember(y_names,'ER40HAP')) = [];
y_types = cell(1, size(y_names, 2));
y_types(:) = {'continuous'};
outname = fullfile(out_dir, 'beh_measures.mat');
delimiter = ',';
y = CBIG_read_y_from_csv(csv_files, subject_header, y_names, y_types, ...
    subject_list, outname, delimiter);
% param.y = y;

subj_list = read_sub_list(subject_list);
ind_nonan = not(any(isnan(y), 2));
subj_list = subj_list(ind_nonan);

% save out subject list after filter behaviors
subject_list_filtered = fullfile(out_dir, 'HCP_diff_roi_subj_list.txt');
fileID = fopen(subject_list_filtered, 'w');
for i = 1:size(subj_list, 2)
    fprintf(fileID, '%s\n', subj_list{i});
end
fclose(fileID);

% save out FC file with filtered subject list
corr_mat = zeros(419, 419, size(subj_list, 2));
for i = 1:size(subj_list, 2)
    ind_fc_file = fullfile(out_dir, 'FC_419',...
        ['FC_419_ROIs_', subj_list{i}, '.mat']);
    tmp = load(ind_fc_file);
    corr_mat(:, :, i) = tmp.corr_mat;
end
save(fullfile(out_dir, 'HCP_FC_S1200_1019.mat'), 'corr_mat');
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
