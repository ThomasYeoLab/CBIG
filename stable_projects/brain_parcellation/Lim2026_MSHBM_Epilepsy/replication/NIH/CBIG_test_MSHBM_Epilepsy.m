% ====================================
%% Help
% ====================================
%
% CBIG_test_MSHBM_Epilepsy
%
% This script tests the MS-HBM model on 14 drug-resistant epilepsy subjects
% from the National Institutes of Health (NIH) dataset using leave-one-run-out
% cross-validation. For each subject, individual-specific parcellations are
% generated from all available training runs (leaving one out at a time) using
% three different group priors: HCP, Du, and NIH (the model trained in Step 1).
% FC homogeneity is computed on the held-out test run for each fold and
% averaged per subject. Group-average parcellations are also evaluated.
% Paired t-tests with BH-FDR correction (q<0.05) compare the six conditions,
% replicating Figure 3A in Lim 2026.
% Last tested on 20 May 2026.
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   Preprocessed NIH fMRI surface files (*sm6.nii.gz) must be present in:
%   $proj_dir/replication/NIH/input/{subject}/surf/
%   Motion censor files must be present in:
%   $proj_dir/replication/NIH/input/{subject}/qc/
%   Group priors (Params_Final.mat) and group parcellation labels (.mat) must
%   be present in $proj_dir/replication/mshbm/{HCP|Du|NIH}/.
%
% OUTPUTS:
%   Homolist pointer files (one per fold per subject):
%   $proj_dir/replication/NIH/results/homo_lists/{subject}/lh_homolist_{bld}.txt
%   Individual parcellations (one per prior per CV fold):
%   $proj_dir/replication/NIH/results/mshbm/{subject}/{prior}_cv_{n}/
%   FC homogeneity results (one per condition per CV fold):
%   $proj_dir/replication/NIH/results/homo/{subject}/{subject}_homo_{condition}_cv_{n}.mat
%   Prints mean +/- SD homogeneity and BH-FDR corrected paired t-test results to console.
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, 'stable_projects/brain_parcellation/Lim2026_MSHBM_Epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');
if isempty(data_dir)
    config_path = fullfile(proj_dir, 'replication/config/CBIG_epilepsy_config.sh');
    [~, extracted] = system(['bash -c "source ' config_path ' && echo $CBIG_EPILEPSY_DATA_DIR"']);
    data_dir = strtrim(extracted);
end
if isempty(data_dir)
    error('CBIG_EPILEPSY_DATA_DIR is not set. Source replication/config/CBIG_epilepsy_config.sh.');
end
analysis_dir  = fullfile(data_dir, 'NIH/input');  % read-only input data
results_dir   = fullfile(proj_dir, 'replication/NIH/results');  % writable outputs
mshbm_ref_dir = fullfile(proj_dir, 'replication/mshbm');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

% Paths to group priors (Params_Final.mat)
prior_path_HCP = fullfile(mshbm_ref_dir, 'HCP/Params_Final.mat');
prior_path_Du  = fullfile(mshbm_ref_dir, 'Du/Params_Final.mat');
prior_path_NIH = fullfile(mshbm_ref_dir, 'NIH/Params_Final.mat');

% Paths to group parcellation labels
grppar_path_HCP = fullfile(mshbm_ref_dir, 'HCP/group.mat');
grppar_path_Du  = fullfile(mshbm_ref_dir, 'Du/group_consensus_raw.mat');
grppar_path_NIH = fullfile(mshbm_ref_dir, 'NIH/group.mat');

% Medial wall mask
medial_wall_path = fullfile(mshbm_ref_dir, 'union_medial_wall.mat');

% Add paths to required functions
addpath(proj_dir);
addpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Load Test Subjects
% ====================================

subject_list_file = fullfile(analysis_dir, 'NIH_test_14_subid.csv');
fid = fopen(subject_list_file);
fgetl(fid);
subject_list_cell = textscan(fid, '%s', 'Delimiter', '\t');
fclose(fid);
subid = subject_list_cell{1};
run_pattern = 'bld([\d]+)';

% ====================================
%% Build Leave-One-Run-Out CV Structure
% ====================================

subj_perm = cell(numel(subid), 3);
for sub_n = 1:numel(subid)
    subject = subid{sub_n};
    surf_files = dir(fullfile(analysis_dir, subject, 'surf', 'lh*sm6.nii.gz'));
    surf_file_names = {};
    for i = 1:numel(surf_files)
        surf_file_names{i} = surf_files(i).name(4:end);
    end

    % Training set for each fold: all runs except the held-out run
    subj_perm{sub_n, 1} = cellfun(@(f) surf_file_names(~strcmp(surf_file_names, f)), ...
        surf_file_names, 'UniformOutput', false);
    % Test run for each fold: the held-out run
    subj_perm{sub_n, 2} = cellfun(@(f) surf_file_names(strcmp(surf_file_names, f)), ...
        surf_file_names, 'UniformOutput', false);
    subj_perm{sub_n, 3} = subject;
end

% Exclude subjects with only 1 run (cannot perform leave-one-run-out CV)
rows_to_remove = cellfun(@(x) numel(x) <= 1, subj_perm(:, 1)) | ...
    cellfun(@(x) numel(x) <= 1, subj_perm(:, 2));
cv_subject = subj_perm(~rows_to_remove, :);
fprintf('CV subjects: %d of %d\n', size(cv_subject, 1), numel(subid));

% ====================================
%% Create Homo Directories and Homolist Files
% ====================================

for sub_n = 1:size(cv_subject, 1)
    subject = cv_subject{sub_n, 3};
    homo_dir = fullfile(results_dir, 'homo_lists', subject);
    if ~exist(homo_dir, 'dir')
        mkdir(homo_dir);
    end

    % Each CV fold: write homolist pointing to the test run fMRI file
    num_folds = numel(cv_subject{sub_n, 2});
    for fold_n = 1:num_folds
        test_file = cv_subject{sub_n, 2}{fold_n}{1};
        bld = regexp(test_file, run_pattern, 'match');
        bld_str = bld{1};

        lh_path = fullfile(analysis_dir, subject, 'surf', ['lh.' test_file]);
        rh_path = fullfile(analysis_dir, subject, 'surf', ['rh.' test_file]);

        lhlist = fullfile(homo_dir, ['lh_homolist_' bld_str '.txt']);
        rhlist = fullfile(homo_dir, ['rh_homolist_' bld_str '.txt']);

        fid = fopen(lhlist, 'w');
        fprintf(fid, '%s\n', lh_path);
        fclose(fid);

        fid = fopen(rhlist, 'w');
        fprintf(fid, '%s\n', rh_path);
        fclose(fid);
    end
end

% ====================================
%% Generate Individual Parcellations (Leave-One-Run-Out CV)
% ====================================

% Three group priors to evaluate
mshbm_flags  = {'HCP', 'Du', 'NIH'};
prior_paths  = {prior_path_HCP, prior_path_Du, prior_path_NIH};

w = '80';
c = '10';

for sub_n = 1:size(cv_subject, 1)
    subject = cv_subject{sub_n, 3};
    num_folds = numel(cv_subject{sub_n, 1});

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};
        prior_path = prior_paths{mshbm_n};

        for fold_n = 1:num_folds
            out_flag = ['_cv_' num2str(fold_n)];
            out_dir = fullfile(results_dir, 'mshbm', subject, [mshbm_flag out_flag]);

            if exist(fullfile(out_dir, 'ind_parcellation/test_set', ...
                    'Ind_parcellation_MSHBM_sub1_w80_MRF10.mat'), 'file')
                fprintf('Parcellation exists: %s %s fold %d, skipping\n', ...
                    subject, mshbm_flag, fold_n);
                continue
            end

            % Training runs for this fold (all runs except the test run)
            train_files = cv_subject{sub_n, 1}{fold_n};
            lh_list = cellfun(@(f) fullfile(analysis_dir, subject, 'surf', ['lh.' f]), ...
                train_files, 'UniformOutput', false);
            rh_list = cellfun(@(f) fullfile(analysis_dir, subject, 'surf', ['rh.' f]), ...
                train_files, 'UniformOutput', false);
            censor_list = cellfun(@(f) get_censor_path(analysis_dir, subject, f, run_pattern), ...
                train_files, 'UniformOutput', false);

            params = struct();
            params.project_dir  = out_dir;
            params.group_prior  = prior_path;
            params.lh_fMRI_list = lh_list;
            params.rh_fMRI_list = rh_list;
            params.censor_list  = censor_list;
            params.target_mesh  = 'fsaverage6';
            params.w            = w;
            params.c            = c;

            CBIG_MSHBM_Epilepsy_LI(params);
            fprintf('Done: %s %s fold %d\n', subject, mshbm_flag, fold_n);
        end
    end
end

% ====================================
%% Compute Homogeneity for Individual Parcellations
% ====================================

load(medial_wall_path, 'mw_lh', 'mw_rh');

for sub_n = 1:size(cv_subject, 1)
    subject = cv_subject{sub_n, 3};
    homo_dir = fullfile(results_dir, 'homo_lists', subject);
    result_homo_dir = fullfile(results_dir, 'homo', subject);
    if ~exist(result_homo_dir, 'dir')
        mkdir(result_homo_dir);
    end
    num_folds = numel(cv_subject{sub_n, 2});

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};

        for fold_n = 1:num_folds
            out_flag = [mshbm_flag '_cv_' num2str(fold_n)];
            homo_file = fullfile(result_homo_dir, [subject '_homo_' out_flag '.mat']);

            if exist(homo_file, 'file')
                fprintf('Homo exists: %s %s, skipping\n', subject, out_flag);
                continue
            end

            parcel_file = fullfile(results_dir, 'mshbm', subject, ...
                [mshbm_flag '_cv_' num2str(fold_n)], ...
                'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

            if ~exist(parcel_file, 'file')
                warning('Parcellation not found: %s', parcel_file);
                continue
            end

            sub_IP = load(parcel_file);
            lh_labels = double(sub_IP.lh_labels);
            rh_labels = double(sub_IP.rh_labels);
            lh_labels(mw_lh) = 0;
            rh_labels(mw_rh) = 0;

            test_file = cv_subject{sub_n, 2}{fold_n}{1};
            bld = regexp(test_file, run_pattern, 'match');
            lhlist = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
            rhlist = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

            homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
                lh_labels, rh_labels, 'fsaverage6', lhlist, rhlist);
            save(homo_file, 'homo_with_weight');
            fprintf('Saved: %s\n', homo_file);
        end
    end
end

% ====================================
%% Compute Homogeneity for Group Parcellations
% ====================================

grppar_paths = {grppar_path_HCP, grppar_path_Du, grppar_path_NIH};
grp_flags    = {'HCP_group', 'Du_group', 'NIH_group'};

for sub_n = 1:size(cv_subject, 1)
    subject = cv_subject{sub_n, 3};
    homo_dir = fullfile(results_dir, 'homo_lists', subject);
    result_homo_dir = fullfile(results_dir, 'homo', subject);
    if ~exist(result_homo_dir, 'dir')
        mkdir(result_homo_dir);
    end
    num_folds = numel(cv_subject{sub_n, 2});

    for parcel_n = 1:numel(grp_flags)
        grp_flag  = grp_flags{parcel_n};
        parcellab = grppar_paths{parcel_n};

        sub_IP = load(parcellab);
        lh_labels_grp = double(sub_IP.lh_labels);
        rh_labels_grp = double(sub_IP.rh_labels);
        lh_labels_grp(mw_lh) = 0;
        rh_labels_grp(mw_rh) = 0;

        for fold_n = 1:num_folds
            out_flag  = [grp_flag '_cv_' num2str(fold_n)];
            homo_file = fullfile(result_homo_dir, [subject '_homo_' out_flag '.mat']);

            if exist(homo_file, 'file')
                fprintf('Homo exists: %s %s, skipping\n', subject, out_flag);
                continue
            end

            test_file = cv_subject{sub_n, 2}{fold_n}{1};
            bld = regexp(test_file, run_pattern, 'match');
            lhlist = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
            rhlist = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

            homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
                lh_labels_grp, rh_labels_grp, 'fsaverage6', lhlist, rhlist);
            save(homo_file, 'homo_with_weight');
            fprintf('Saved: %s\n', homo_file);
        end
    end
end

% ====================================
%% Load Homogeneity and Aggregate Per Subject
% ====================================

% Six conditions: HCP group, HCP individual, Du group, Du individual, NIH group, NIH individual
parcel_patterns = {'HCP_group_cv_*', 'HCP_cv_*', ...
    'Du_group_cv_*', 'Du_cv_*', ...
    'NIH_group_cv_*', 'NIH_cv_*'};
parcel_names = {'HCP_group', 'HCP_IP', ...
    'Du_group', 'Du_IP', ...
    'NIH_group', 'NIH_IP'};

sub_homo = NaN(size(cv_subject, 1), numel(parcel_names));

for parcel_n = 1:numel(parcel_patterns)
    for sub_n = 1:size(cv_subject, 1)
        subject = cv_subject{sub_n, 3};
        result_homo_dir = fullfile(results_dir, 'homo', subject);
        homo_files = dir(fullfile(result_homo_dir, ...
            ['*_homo_' parcel_patterns{parcel_n} '.mat']));

        if numel(homo_files) == 0
            warning('No homo files for %s %s', subject, parcel_names{parcel_n});
            continue
        end

        fold_homos = NaN(numel(homo_files), 1);
        for f = 1:numel(homo_files)
            tmp = load(fullfile(homo_files(f).folder, homo_files(f).name));
            fold_homos(f) = tmp.homo_with_weight;
        end
        sub_homo(sub_n, parcel_n) = mean(fold_homos, 'omitnan');
    end
end

% ====================================
%% Paired T-Tests and FDR Correction
% ====================================

p_vals = [];
pair_names = {};
for i = 1:numel(parcel_names)
    for j = (i+1):numel(parcel_names)
        [~, p] = ttest(sub_homo(:, i), sub_homo(:, j));
        p_vals(end+1)     = p;
        pair_names{end+1} = [parcel_names{i} ' vs ' parcel_names{j}];
    end
end

p_fdr = mafdr(p_vals, 'BHFDR', true);

% ====================================
%% Print Results
% ====================================

fprintf('\n--- Mean FC Homogeneity per Condition ---\n');
for parcel_n = 1:numel(parcel_names)
    fprintf('%s: %.4f +/- %.4f\n', parcel_names{parcel_n}, ...
        nanmean(sub_homo(:, parcel_n)), nanstd(sub_homo(:, parcel_n)));
end

fprintf('\n--- Paired T-Test Results (FDR corrected, q<0.05) ---\n');
for k = 1:numel(pair_names)
    sig = '';
    if p_fdr(k) < 0.05
        sig = ' *';
    end
    fprintf('%s: p_fdr=%.4f%s\n', pair_names{k}, p_fdr(k), sig);
end

% ====================================
%% Cleanup Paths
% ====================================

rmpath(proj_dir);
rmpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Local Functions
% ====================================

function censor_path = get_censor_path(analysis_dir, subject, surf_file, run_pattern)
% Returns the censor file path corresponding to a given surface fMRI filename.
bld = regexp(surf_file, run_pattern, 'match');
censor_path = fullfile(analysis_dir, subject, 'qc', ...
    [subject '_' bld{1} '_FDRMS0.2_DVARS50_motion_outliers.txt']);
end
