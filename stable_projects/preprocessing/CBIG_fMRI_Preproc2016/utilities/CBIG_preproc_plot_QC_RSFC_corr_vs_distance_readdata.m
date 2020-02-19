function [QC, RSFC, y_label] = CBIG_preproc_plot_QC_RSFC_corr_vs_distance_readdata( ...
    sub_dir, sub_list, QC_type, RSFC_stem )

% CBIG_preproc_plot_QC_RSFC_corr_vs_distance_readdata( sub_dir, sub_list, QC_type, RSFC_stem )
% 
% Read QC (quality control) measure, RSFC (resting-state functional
% connectivity) matrix of each subject, and return a QC vector of all
% subjects and a 3D RSFC matrix of all subjects.
% 
% This code assumes the data went through our CBIG_fMRI_Preproc2016
% preprocessing pipeline. For example, the functional connecitivity matrix
% for each subject is
% <path_to_subject>/<subject_name>/FC_metrics/Pearson_r/<subjact_name>_<RSFC_stem>.mat;
% the mean FD for each subject is
% <path_to_subject>/<subject_name>/bold/<run#>/<subject_name>_bld<run#>_rest_skip4_stc_mc_rel_mean.rms;
% etc.
% 
% Inputs:
%     - sub_dir:
%       The absolute path where the preprocessed subjects are stored.
% 
%     - sub_list:
%       The full name of a text file, in which each line is one subject ID.
% 
%     - QC_type:
%       A string specify what is the QC (quality control) measure. Currenly
%       we only support for 'FD_mean', 'FD_std', and 'censored_frames'.
%       'FD_mean' means the QC measure is the mean framewise displacement.
%       'FD_std' means the QC measure is the standard deviation of
%       framewise displacement. 'censored_frames' means the number of
%       censored frames determined by FD and/or DVARS.
% 
%     - RSFC_stem:
%       The stem of RSFC (resting-state functional connectivity) filenames.
%       For example, if the RSFC filename for one subject is:
%       <path_to_data>/Sub0001_Ses1_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_all2all.mat,
%       then RSFC_stem =
%       '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_all2all'.
% 
% Outputs:
%     - QC:
%       An N x 1 vector contains QC measures of all subjects (e.g. mean
%       framewise displacement), where N is number of subjects.
% 
%     - RSFC:
%       An M x M x N matrix of functional connectivity for all subjects,
%       where M is number of ROIs.
% 
%     - y_label:
%       The ylabel that can be displayed on the QC-FC plot, if the users
%       will call CBIG_preproc_plot_QC_RSFC_corr_vs_distance_matrix.m
%       later. For example: 'FD (mean) - RSFC correlation'.
% 
% Example 1:
%     [QC, RSFC, y_label] = CBIG_preproc_plot_QC_RSFC_corr_vs_distance_readdata( ...
% fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', 'preprocessing', 'CBIG_fMRI_Preproc2016', ...
% '100subjects_clustering', 'preproc_out'), fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
% 'CBIG_fMRI_Preproc2016', 'unit_tests', '100subjects_clustering', 'GSP_80_low_motion+20_w_censor.txt'), ...
% 'FD_mean', '_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_all2all' )
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% Parse subject names
fid = fopen(sub_list);
subjects = textscan(fid, '%s\n');
subjects = subjects{1};
num_sub = length(subjects);
fclose(fid);


%% Read in RSFC
% It assumes the user has run CBIG_preproc_FC_metrics.csh to generate RSFC
% file and put into <sub_dir>/<subject_name>/FC_metrics/Pearson_r folder.
RSFC = [];
for i = 1:num_sub
    RSFC_file = fullfile(sub_dir, subjects{i}, 'FC_metrics', 'Pearson_r', [subjects{i} RSFC_stem '.mat']);
    load(RSFC_file)
    RSFC = cat(3, RSFC, corr_mat);
end


%% Different types of QC measures 
QC = zeros(num_sub, 1);
for i = 1:num_sub
    % read run numbers
    boldfile = fullfile(sub_dir, subjects{i}, 'logs', [subjects{i} '.bold']);
    fid = fopen(boldfile, 'r');
    xdata = fgets(fid);
    fclose(fid);
    C = strsplit(xdata, ' ');                % split this line by spaces
    C{end} = strrep(C{end}, char(10), '');   % remove "line feed" charactor
    runs = C(~cellfun('isempty', C));        % remove empty cell arrays
    
    % read qc for each run
    qc = zeros(length(runs), 1);
    
    if(strcmp(QC_type, 'FD_mean'))
        for j = 1:length(runs)
            QC_file = fullfile(sub_dir, subjects{i}, 'bold', runs{j}, [subjects{i} '_bld' runs{j} ...
                '_rest_skip4_stc_mc_rel_mean.rms']);
            qc(j) = dlmread(QC_file);
        end
        QC(i) = mean(qc);
        y_label = 'FD (mean) - RSFC correlation';
        
    elseif(strcmp(QC_type, 'FD_std'))
        for j = 1:length(runs)
            QC_file = fullfile(sub_dir , subjects{i}, 'bold', runs{j}, [subjects{i} '_bld' runs{j} ...
                '_rest_skip4_stc_mc_rel.rms']);
            qc(j) = std(dlmread(QC_file));
        end
        QC(i) = mean(qc);
        y_label = 'FD (std) - RSFC correlation';
        
    elseif(strcmp(QC_type, 'censored_frames'))
        for j = 1:length(runs)
            QC_file = fullfile(sub_dir, subjects{i}, 'qc', [subjects{i} '_bld' runs{j} ...
                '_FDRMS0.2_DVARS50_motion_outliers.txt']);
            outliers = dlmread(QC_file);
            qc(j) = sum(outliers==0);
        end
        QC(i) = sum(qc) / length(runs);
        y_label = '# Censored frames - RSFC correlation';
        
    else
        error('Unknown QC measure type.')
    end
end



end
