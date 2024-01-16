function CBIG_MM_HCP_generate_FC(in_dir, out_dir)

% CBIG_MM_HCP_generate_FC(in_dir, out_dir)
% 
% This function generate functional connectivity data for HCP S1200
%
% Inputs:
%   - in_dir
%     Path of the input time series data for HCP S1200.
%
%   - out_dir
%     Path of the your output functional connectivity data for HCP
%     S1200.
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if ~exist(out_dir, 'dir')
   mkdir(out_dir)
end

%% Read info in input files
% read subject run list
dir_inf = struct2cell(dir(fullfile(in_dir, '*.mat')));
files = dir_inf(1, :);
subject_list = unique(cellfun(@(x) str2num(x(13:18)), files, 'UniformOutput', true));
disp(['we have ', num2str(length(subject_list)), ' subjects with time series'])

for i = 1:length(subject_list)
    subject = subject_list(i);
    subject_files = struct2cell(dir(fullfile(in_dir, ['*', num2str(subject), '*.mat'])));
    subject_folder = subject_files(2, :);
    subject_files = subject_files(1, :);
    for j = 1:length(subject_files)
        file = [fullfile(subject_folder{j}, subject_files{j})];
        tc = load(file);
        tc = tc.TC;
        tc = bsxfun(@minus, tc, mean(tc, 1));
        tc = bsxfun(@times, tc, 1./sqrt(sum(tc.^2, 1)));
        subj_corr_mat = tc' * tc;
        if j == 1
            subj_z_mat = CBIG_StableAtanh(subj_corr_mat); % Fisher r-to-z transform
        else
            subj_z_mat = subj_z_mat + CBIG_StableAtanh(subj_corr_mat);
        end
    end
    subj_z_mat = subj_z_mat/length(subject_files); % average across number of runs
    corr_mat = tanh(subj_z_mat);
    save_name = fullfile(out_dir, ['/FC_419_ROIs_', num2str(subject), '.mat']);
    save(save_name,'corr_mat');
end
end
