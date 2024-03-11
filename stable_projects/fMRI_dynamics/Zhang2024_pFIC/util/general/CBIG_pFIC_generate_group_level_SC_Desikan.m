function group_level_SC = CBIG_pFIC_generate_group_level_SC_Desikan(subject_list, subject_level_SC_path)

% group_level_SC = CBIG_pFIC_generate_group_level_SC_Desikan(subject_list, subject_level_SC_path)
% This function generates a group-averaged structural connectivitiy (SC),
% the subject-level SC is assumed to be computed using MRtrix iFOD2
% tractography with SIFT2 filtering. Entries that are absent from more 
% than 50% of the subjects are masked out. 
% Input:
%   - subject_list: absolute path to an n-by-1 vector where each row is a subject ID (n = number of subjects)
%   - subject_level_SC_path: absolute path to the directory with subject-level SC
% Output:
%   - group_level_SC: a 68-by-68 group-averaged SC. Here, we assume using the Desikan parcellation.
% Example:
% SC_train = CBIG_pFIC_generate_group_level_SC_Desikan('Zhang2024_pFIC/replication/HCP/input/training_subject_list.txt',
%   'tractography/iFOD2/');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% read subject list
fid = fopen(subject_list, 'r');
stringPattern = textscan(fid,'%s');
fclose(fid);
subject_list = string(stringPattern{:});

group_level_SC = zeros(68);
SC_subject_all = zeros(68, 68, length(subject_list));

for i = 1:length(subject_list)
   subject = subject_list{i}; 
   SC_subject = csvread([subject_level_SC_path '/' subject '/connectomes/connectome_DK_82Parcels_SIFT2.csv']);
   SC_subject(35:48, :) = []; % exlcuding subcortical, refer to /Desikan_82_SC.txt for more details
   SC_subject(:, 35:48) = [];
   SC_subject(SC_subject < 0) = 0;
   SC_subject_all(:, :, i) = SC_subject; 
end

% keep track of how many subjects have zero in an entry
zero_count_map = sum(SC_subject_all == 0, 3);

% mask out entries that are absent from more than 50% of the subjects 
SC_mask = zero_count_map <= length(subject_list)/2;
SC_subject_all = SC_mask .* SC_subject_all;

% take natural log of each entry and average the non-zeros values across
% subjects
for i = 1:68
    for j = 1:68
        SC_entry_all = SC_subject_all(i, j, :);
        SC_entry_all = SC_entry_all(:);
        SC_entry_all_log = log(SC_entry_all);
        SC_entry_all_log(isinf(SC_entry_all_log)) = 0;
        group_level_SC(i, j) = mean(nonzeros(SC_entry_all_log));
    end
end
group_level_SC(isnan(group_level_SC)) = 0;

% the main diagonal is set to 0 (i.e. no recurrent connections)
group_level_SC(logical(eye(68))) = 0;

end
