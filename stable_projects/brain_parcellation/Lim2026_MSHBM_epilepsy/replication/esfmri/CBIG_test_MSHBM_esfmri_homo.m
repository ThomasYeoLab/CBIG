% ====================================
%% Help
% ====================================
%
% CBIG_test_MSHBM_esfmri_homo
%
% This script tests the MS-HBM model on the electrostimulation fMRI
% (esfMRI) dataset using leave-one-run-out cross-validation. Two
% independent homogeneity tests are performed:
%
%   (1) Pre-op resting state (rs_pp): 20 subjects, surf files are
%       motion-censored resting-state scans (*fs6_sm6_censored.nii.gz).
%       Only runs passing mc_test (FDRMS0.5_DVARS50; >= 2 valid runs per
%       subject) are included per run.
%
%   (2) Post-op no-stim periods (es_pp): 17 subjects, surf files are
%       pre-filtered no-stimulation epochs (*_nostim_discard.nii.gz).
%       mc_test is used only as a subject-level inclusion threshold
%       (>= 2 valid runs); all nostim_discard files for included subjects
%       are used. censor_list is set to NONE because files are
%       already frame-selected.
%
% For each dataset, individual-specific parcellations are generated from
% training runs (leaving one out at a time) using three group priors: HCP,
% Du, and NIH (NIH_34). FC homogeneity is computed on the held-out
% run and averaged per subject. Group-average parcellations are also
% evaluated. Paired t-tests with BH-FDR correction (q<0.05) compare the
% six conditions, replicating the esfMRI validation results in Lim 2026.
% Last tested on 20 May 2026.
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   Preprocessed resting-state surface files (*fs6_sm6_censored.nii.gz) must be in:
%   $proj_dir/replication/esfmri/input/rs_pp/{subject}/surf/
%   Preprocessed no-stim surface files (*_nostim_discard.nii.gz) must be in:
%   $proj_dir/replication/esfmri/input/es_pp/{subject}/surf/
%   Motion QC files (mc_test.mat) must be present in the rs_pp and es_pp directories.
%   Group priors (Params_Final.mat) and group parcellation labels (.mat) must
%   be present in $proj_dir/replication/mshbm/{HCP|Du|NIH}/.
%
% OUTPUTS:
%   Homolist pointer files (one per fold per subject):
%   $proj_dir/replication/esfmri/results/homo_lists/{rs_pp|es_pp}/{subject}/lh_homolist_{bld}.txt
%   Individual parcellations (one per prior per CV fold, for rs_pp and es_pp):
%   $proj_dir/replication/esfmri/results/mshbm/{rs_pp|es_pp}/{subject}/{prior}_cv_{n}/
%   FC homogeneity results (one per condition per CV fold):
%   $proj_dir/replication/esfmri/results/homo/{rs_pp|es_pp}/{subject}/
%       {subject}_homo_{condition}_cv_{n}.mat
%   Prints mean +/- SD homogeneity and BH-FDR corrected paired t-test results to console.
%   $proj_dir/replication/esfmri/results/homo_rs_pp_summary.csv  (mean +/- SD per condition)
%   $proj_dir/replication/esfmri/results/homo_rs_pp_stats.csv    (pairwise t-test results)
%   $proj_dir/replication/esfmri/results/homo_es_pp_summary.csv  (mean +/- SD per condition)
%   $proj_dir/replication/esfmri/results/homo_es_pp_stats.csv    (pairwise t-test results)
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, 'stable_projects/brain_parcellation/Lim2026_MSHBM_epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');
if isempty(data_dir)
    config_path = fullfile(proj_dir, 'replication/config/CBIG_epilepsy_config.sh');
    [~, extracted] = system(['bash -c "source ' config_path ' && echo $CBIG_EPILEPSY_DATA_DIR"']);
    data_dir = strtrim(extracted);
end
if isempty(data_dir)
    error('CBIG_EPILEPSY_DATA_DIR is not set. Source replication/config/CBIG_epilepsy_config.sh.');
end
input_dir     = fullfile(data_dir, 'esfmri/input');  % read-only input data
results_dir   = fullfile(proj_dir, 'replication/esfmri/results');  % writable outputs
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
load(medial_wall_path, 'mw_lh', 'mw_rh');

% Three group priors and corresponding group parcellation labels
mshbm_flags  = {'HCP', 'Du', 'NIH'};
prior_paths  = {prior_path_HCP, prior_path_Du, prior_path_NIH};
grppar_paths = {grppar_path_HCP, grppar_path_Du, grppar_path_NIH};
grp_flags    = {'HCP_group', 'Du_group', 'NIH_group'};

% Six condition labels (group and individual per prior)
parcel_patterns = {'HCP_group_cv_*', 'HCP_cv_*', ...
    'Du_group_cv_*', 'Du_cv_*', ...
    'NIH_group_cv_*', 'NIH_cv_*'};
parcel_names = {'HCP_group', 'HCP_IP', ...
    'Du_group', 'Du_IP', ...
    'NIH_group', 'NIH_IP'};

w = '80';
c = '10';
run_pattern = 'bld([\d]+)';

% Add paths to required functions
addpath(proj_dir);
addpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% ===== Part 1: Pre-op Resting State (rs_pp) =====
% ====================================

fprintf('\n====== Part 1: rs_pp (pre-op resting state) ======\n');

% ====================================
%% rs_pp: Load Subjects and Motion QC
% ====================================

rspp_dir  = fullfile(input_dir, 'rs_pp');
subid_file = fullfile(input_dir, 'rsfmri_subid.csv');
fid = fopen(subid_file);
fgetl(fid);
subid_rs = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid_rs = subid_rs{1};

mc_rs = load(fullfile(rspp_dir, 'mc_test.mat'), 'mc_test');
mc_test_rs = mc_rs.mc_test;

% ====================================
%% rs_pp: Build Leave-One-Run-Out CV Structure
% ====================================

% Only runs passing mc_test (pass == 1) are included per subject.
% Subjects with fewer than 2 valid runs are excluded.
subj_perm_rs = cell(numel(subid_rs), 3);
for sub_n = 1:numel(subid_rs)
    subject = subid_rs{sub_n};
    subj    = strrep(subject, '-', '_');

    if ~isfield(mc_test_rs, subj) || sum(mc_test_rs.(subj).pass) < 2
        continue
    end

    surf_files = dir(fullfile(rspp_dir, subject, 'surf', ...
        'lh.*fs6_sm6_censored.nii.gz'));
    surf_file_names = {};
    for i = 1:numel(surf_files)
        if mc_test_rs.(subj).pass(i) == 0
            continue
        end
        surf_file_names{end+1} = surf_files(i).name(4:end); %#ok<SAGROW>
    end

    subj_perm_rs{sub_n, 1} = cellfun(@(f) surf_file_names(~strcmp(surf_file_names, f)), ...
        surf_file_names, 'UniformOutput', false);
    subj_perm_rs{sub_n, 2} = cellfun(@(f) surf_file_names(strcmp(surf_file_names, f)), ...
        surf_file_names, 'UniformOutput', false);
    subj_perm_rs{sub_n, 3} = subject;
end

rows_to_remove = cellfun(@(x) numel(x) <= 1, subj_perm_rs(:, 1)) | ...
    cellfun(@(x) numel(x) <= 1, subj_perm_rs(:, 2));
cv_rs = subj_perm_rs(~rows_to_remove, :);
fprintf('rs_pp CV subjects: %d of %d\n', size(cv_rs, 1), numel(subid_rs));

% ====================================
%% rs_pp: Create Homo Directories and Homolist Files
% ====================================

for sub_n = 1:size(cv_rs, 1)
    subject  = cv_rs{sub_n, 3};
    homo_dir = fullfile(results_dir, 'homo_lists', 'rs_pp', subject);
    if ~exist(homo_dir, 'dir')
        mkdir(homo_dir);
    end

    num_folds = numel(cv_rs{sub_n, 2});
    for fold_n = 1:num_folds
        test_file = cv_rs{sub_n, 2}{fold_n}{1};
        bld       = regexp(test_file, run_pattern, 'match');
        bld_str   = bld{1};

        lh_path = fullfile(rspp_dir, subject, 'surf', ['lh.' test_file]);
        rh_path = fullfile(rspp_dir, subject, 'surf', ['rh.' test_file]);

        fid = fopen(fullfile(homo_dir, ['lh_homolist_' bld_str '.txt']), 'w');
        fprintf(fid, '%s\n', lh_path);
        fclose(fid);

        fid = fopen(fullfile(homo_dir, ['rh_homolist_' bld_str '.txt']), 'w');
        fprintf(fid, '%s\n', rh_path);
        fclose(fid);
    end
end

% ====================================
%% rs_pp: Generate Individual Parcellations
% ====================================

for sub_n = 1:size(cv_rs, 1)
    subject   = cv_rs{sub_n, 3};
    num_folds = numel(cv_rs{sub_n, 1});

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};
        prior_path = prior_paths{mshbm_n};

        for fold_n = 1:num_folds
            out_dir = fullfile(results_dir, 'mshbm', 'rs_pp', subject, ...
                [mshbm_flag '_cv_' num2str(fold_n)]);

            if exist(fullfile(out_dir, 'ind_parcellation/test_set', ...
                    'Ind_parcellation_MSHBM_sub1_w80_MRF10.mat'), 'file')
                fprintf('Exists: %s %s fold %d, skipping\n', subject, mshbm_flag, fold_n);
                continue
            end

            train_files = cv_rs{sub_n, 1}{fold_n};
            lh_list     = cellfun(@(f) fullfile(rspp_dir, subject, 'surf', ['lh.' f]), ...
                train_files, 'UniformOutput', false);
            rh_list     = cellfun(@(f) fullfile(rspp_dir, subject, 'surf', ['rh.' f]), ...
                train_files, 'UniformOutput', false);

            params              = struct();
            params.project_dir  = out_dir;
            params.group_prior  = prior_path;
            params.lh_fMRI_list = lh_list;
            params.rh_fMRI_list = rh_list;
            params.censor_list  = 'NONE';
            params.target_mesh  = 'fsaverage6';
            params.w            = w;
            params.c            = c;

            CBIG_MSHBM_Epilepsy_LI(params);
            fprintf('Done: %s %s fold %d\n', subject, mshbm_flag, fold_n);
        end
    end
end

% ====================================
%% rs_pp: Compute Homogeneity for Individual Parcellations
% ====================================

for sub_n = 1:size(cv_rs, 1)
    subject         = cv_rs{sub_n, 3};
    homo_dir        = fullfile(results_dir, 'homo_lists', 'rs_pp', subject);
    result_homo_dir = fullfile(results_dir, 'homo', 'rs_pp', subject);
    if ~exist(result_homo_dir, 'dir')
        mkdir(result_homo_dir);
    end
    num_folds = numel(cv_rs{sub_n, 2});

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};

        for fold_n = 1:num_folds
            out_flag  = [mshbm_flag '_cv_' num2str(fold_n)];
            homo_file = fullfile(result_homo_dir, [subject '_homo_' out_flag '.mat']);

            if exist(homo_file, 'file')
                fprintf('Homo exists: %s %s, skipping\n', subject, out_flag);
                continue
            end

            parcel_file = fullfile(results_dir, 'mshbm', 'rs_pp', subject, ...
                [mshbm_flag '_cv_' num2str(fold_n)], ...
                'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

            if ~exist(parcel_file, 'file')
                warning('Parcellation not found: %s', parcel_file);
                continue
            end

            sub_IP    = load(parcel_file);
            lh_labels = double(sub_IP.lh_labels);
            rh_labels = double(sub_IP.rh_labels);
            lh_labels(mw_lh) = 0;
            rh_labels(mw_rh) = 0;

            test_file = cv_rs{sub_n, 2}{fold_n}{1};
            bld       = regexp(test_file, run_pattern, 'match');
            lhlist    = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
            rhlist    = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

            homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
                lh_labels, rh_labels, 'fsaverage6', lhlist, rhlist);
            save(homo_file, 'homo_with_weight');
            fprintf('Saved: %s\n', homo_file);
        end
    end
end

% ====================================
%% rs_pp: Compute Homogeneity for Group Parcellations
% ====================================

for sub_n = 1:size(cv_rs, 1)
    subject         = cv_rs{sub_n, 3};
    homo_dir        = fullfile(results_dir, 'homo_lists', 'rs_pp', subject);
    result_homo_dir = fullfile(results_dir, 'homo', 'rs_pp', subject);
    if ~exist(result_homo_dir, 'dir')
        mkdir(result_homo_dir);
    end
    num_folds = numel(cv_rs{sub_n, 2});

    for parcel_n = 1:numel(grp_flags)
        grp_flag  = grp_flags{parcel_n};
        sub_IP        = load(grppar_paths{parcel_n});
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

            test_file = cv_rs{sub_n, 2}{fold_n}{1};
            bld       = regexp(test_file, run_pattern, 'match');
            lhlist    = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
            rhlist    = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

            homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
                lh_labels_grp, rh_labels_grp, 'fsaverage6', lhlist, rhlist);
            save(homo_file, 'homo_with_weight');
            fprintf('Saved: %s\n', homo_file);
        end
    end
end

% ====================================
%% rs_pp: Aggregate, Statistics, and Print Results
% ====================================

sub_homo_rs = NaN(size(cv_rs, 1), numel(parcel_names));
for parcel_n = 1:numel(parcel_patterns)
    for sub_n = 1:size(cv_rs, 1)
        subject         = cv_rs{sub_n, 3};
        result_homo_dir = fullfile(results_dir, 'homo', 'rs_pp', subject);
        homo_files      = dir(fullfile(result_homo_dir, ...
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
        sub_homo_rs(sub_n, parcel_n) = mean(fold_homos, 'omitnan');
    end
end

p_vals_rs = [];
pair_names_rs = {};
for i = 1:numel(parcel_names)
    for j = (i+1):numel(parcel_names)
        [~, p]               = ttest(sub_homo_rs(:, i), sub_homo_rs(:, j));
        p_vals_rs(end+1)     = p; %#ok<SAGROW>
        pair_names_rs{end+1} = [parcel_names{i} ' vs ' parcel_names{j}]; %#ok<SAGROW>
    end
end
p_fdr_rs = mafdr(p_vals_rs, 'BHFDR', true);

fprintf('\n--- rs_pp Mean FC Homogeneity per Condition ---\n');
for parcel_n = 1:numel(parcel_names)
    fprintf('%s: %.4f +/- %.4f\n', parcel_names{parcel_n}, ...
        nanmean(sub_homo_rs(:, parcel_n)), nanstd(sub_homo_rs(:, parcel_n)));
end

fprintf('\n--- rs_pp Paired T-Test Results (FDR corrected, q<0.05) ---\n');
for k = 1:numel(pair_names_rs)
    sig = '';
    if p_fdr_rs(k) < 0.05
        sig = ' *';
    end
    fprintf('%s: p_fdr=%.4f%s\n', pair_names_rs{k}, p_fdr_rs(k), sig);
end

% Save rs_pp results to CSV
summary_file = fullfile(results_dir, 'homo_rs_pp_summary.csv');
fid = fopen(summary_file, 'w');
fprintf(fid, 'condition,mean_homo,sd_homo\n');
for p = 1:numel(parcel_names)
    fprintf(fid, '%s,%.6f,%.6f\n', parcel_names{p}, ...
        nanmean(sub_homo_rs(:, p)), nanstd(sub_homo_rs(:, p)));
end
fclose(fid);
fprintf('\nSaved: %s\n', summary_file);

stats_file = fullfile(results_dir, 'homo_rs_pp_stats.csv');
fid = fopen(stats_file, 'w');
fprintf(fid, 'pair,p_uncorr,p_fdr,sig_uncorr,sig_fdr\n');
for k = 1:numel(pair_names_rs)
    fprintf(fid, '%s,%.6f,%.6f,%d,%d\n', pair_names_rs{k}, ...
        p_vals_rs(k), p_fdr_rs(k), p_vals_rs(k) < 0.05, p_fdr_rs(k) < 0.05);
end
fclose(fid);
fprintf('Saved: %s\n', stats_file);

% ====================================
%% ===== Part 2: Post-op No-Stim Periods (es_pp) =====
% ====================================

fprintf('\n====== Part 2: es_pp (post-op no-stim periods) ======\n');

% ====================================
%% es_pp: Load Subjects and Motion QC
% ====================================

espp_dir   = fullfile(input_dir, 'es_pp');
subid_file = fullfile(input_dir, 'esfmri_subid.csv');
fid = fopen(subid_file);
fgetl(fid);
subid_es = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid_es = subid_es{1};

mc_es = load(fullfile(espp_dir, 'mc_test.mat'), 'mc_test');
mc_test_es = mc_es.mc_test;

% ====================================
%% es_pp: Build Leave-One-Run-Out CV Structure
% ====================================

% mc_test is used only as a subject-level threshold (>= 2 valid runs).
% All nostim_discard files for qualifying subjects are included in CV —
% no per-run mc_test filtering — because the files are already frame-selected.
subj_perm_es = cell(numel(subid_es), 3);
for sub_n = 1:numel(subid_es)
    subject = subid_es{sub_n};
    subj    = strrep(subject, '-', '_');

    if ~isfield(mc_test_es, subj) || sum(mc_test_es.(subj).pass) < 2
        continue
    end

    surf_files = dir(fullfile(espp_dir, subject, 'surf', ...
        'lh.*_nostim_discard.nii.gz'));
    surf_file_names = {};
    for i = 1:numel(surf_files)
        surf_file_names{end+1} = surf_files(i).name(4:end); %#ok<SAGROW>
    end

    subj_perm_es{sub_n, 1} = cellfun(@(f) surf_file_names(~strcmp(surf_file_names, f)), ...
        surf_file_names, 'UniformOutput', false);
    subj_perm_es{sub_n, 2} = cellfun(@(f) surf_file_names(strcmp(surf_file_names, f)), ...
        surf_file_names, 'UniformOutput', false);
    subj_perm_es{sub_n, 3} = subject;
end

rows_to_remove = cellfun(@(x) numel(x) <= 1, subj_perm_es(:, 1)) | ...
    cellfun(@(x) numel(x) <= 1, subj_perm_es(:, 2));
cv_es = subj_perm_es(~rows_to_remove, :);
fprintf('es_pp CV subjects: %d of %d\n', size(cv_es, 1), numel(subid_es));

% ====================================
%% es_pp: Create Homo Directories and Homolist Files
% ====================================

for sub_n = 1:size(cv_es, 1)
    subject  = cv_es{sub_n, 3};
    homo_dir = fullfile(results_dir, 'homo_lists', 'es_pp', subject);
    if ~exist(homo_dir, 'dir')
        mkdir(homo_dir);
    end

    num_folds = numel(cv_es{sub_n, 2});
    for fold_n = 1:num_folds
        test_file = cv_es{sub_n, 2}{fold_n}{1};
        bld       = regexp(test_file, run_pattern, 'match');
        bld_str   = bld{1};

        lh_path = fullfile(espp_dir, subject, 'surf', ['lh.' test_file]);
        rh_path = fullfile(espp_dir, subject, 'surf', ['rh.' test_file]);

        fid = fopen(fullfile(homo_dir, ['lh_homolist_' bld_str '.txt']), 'w');
        fprintf(fid, '%s\n', lh_path);
        fclose(fid);

        fid = fopen(fullfile(homo_dir, ['rh_homolist_' bld_str '.txt']), 'w');
        fprintf(fid, '%s\n', rh_path);
        fclose(fid);
    end
end

% ====================================
%% es_pp: Generate Individual Parcellations
% ====================================

% censor_list is NONE because nostim_discard files are already frame-selected.
for sub_n = 1:size(cv_es, 1)
    subject   = cv_es{sub_n, 3};
    num_folds = numel(cv_es{sub_n, 1});

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};
        prior_path = prior_paths{mshbm_n};

        for fold_n = 1:num_folds
            out_dir = fullfile(results_dir, 'mshbm', 'es_pp', subject, ...
                [mshbm_flag '_cv_' num2str(fold_n)]);

            if exist(fullfile(out_dir, 'ind_parcellation/test_set', ...
                    'Ind_parcellation_MSHBM_sub1_w80_MRF10.mat'), 'file')
                fprintf('Exists: %s %s fold %d, skipping\n', subject, mshbm_flag, fold_n);
                continue
            end

            train_files = cv_es{sub_n, 1}{fold_n};
            lh_list     = cellfun(@(f) fullfile(espp_dir, subject, 'surf', ['lh.' f]), ...
                train_files, 'UniformOutput', false);
            rh_list     = cellfun(@(f) fullfile(espp_dir, subject, 'surf', ['rh.' f]), ...
                train_files, 'UniformOutput', false);

            params              = struct();
            params.project_dir  = out_dir;
            params.group_prior  = prior_path;
            params.lh_fMRI_list = lh_list;
            params.rh_fMRI_list = rh_list;
            params.censor_list  = 'NONE';
            params.target_mesh  = 'fsaverage6';
            params.w            = w;
            params.c            = c;

            CBIG_MSHBM_Epilepsy_LI(params);
            fprintf('Done: %s %s fold %d\n', subject, mshbm_flag, fold_n);
        end
    end
end

% ====================================
%% es_pp: Compute Homogeneity for Individual Parcellations
% ====================================

for sub_n = 1:size(cv_es, 1)
    subject         = cv_es{sub_n, 3};
    homo_dir        = fullfile(results_dir, 'homo_lists', 'es_pp', subject);
    result_homo_dir = fullfile(results_dir, 'homo', 'es_pp', subject);
    if ~exist(result_homo_dir, 'dir')
        mkdir(result_homo_dir);
    end
    num_folds = numel(cv_es{sub_n, 2});

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};

        for fold_n = 1:num_folds
            out_flag  = [mshbm_flag '_cv_' num2str(fold_n)];
            homo_file = fullfile(result_homo_dir, [subject '_homo_' out_flag '.mat']);

            if exist(homo_file, 'file')
                fprintf('Homo exists: %s %s, skipping\n', subject, out_flag);
                continue
            end

            parcel_file = fullfile(results_dir, 'mshbm', 'es_pp', subject, ...
                [mshbm_flag '_cv_' num2str(fold_n)], ...
                'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

            if ~exist(parcel_file, 'file')
                warning('Parcellation not found: %s', parcel_file);
                continue
            end

            sub_IP    = load(parcel_file);
            lh_labels = double(sub_IP.lh_labels);
            rh_labels = double(sub_IP.rh_labels);
            lh_labels(mw_lh) = 0;
            rh_labels(mw_rh) = 0;

            test_file = cv_es{sub_n, 2}{fold_n}{1};
            bld       = regexp(test_file, run_pattern, 'match');
            lhlist    = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
            rhlist    = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

            homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
                lh_labels, rh_labels, 'fsaverage6', lhlist, rhlist);
            save(homo_file, 'homo_with_weight');
            fprintf('Saved: %s\n', homo_file);
        end
    end
end

% ====================================
%% es_pp: Compute Homogeneity for Group Parcellations
% ====================================

for sub_n = 1:size(cv_es, 1)
    subject         = cv_es{sub_n, 3};
    homo_dir        = fullfile(results_dir, 'homo_lists', 'es_pp', subject);
    result_homo_dir = fullfile(results_dir, 'homo', 'es_pp', subject);
    if ~exist(result_homo_dir, 'dir')
        mkdir(result_homo_dir);
    end
    num_folds = numel(cv_es{sub_n, 2});

    for parcel_n = 1:numel(grp_flags)
        grp_flag  = grp_flags{parcel_n};
        sub_IP        = load(grppar_paths{parcel_n});
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

            test_file = cv_es{sub_n, 2}{fold_n}{1};
            bld       = regexp(test_file, run_pattern, 'match');
            lhlist    = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
            rhlist    = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

            homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
                lh_labels_grp, rh_labels_grp, 'fsaverage6', lhlist, rhlist);
            save(homo_file, 'homo_with_weight');
            fprintf('Saved: %s\n', homo_file);
        end
    end
end

% ====================================
%% es_pp: Aggregate, Statistics, and Print Results
% ====================================

sub_homo_es = NaN(size(cv_es, 1), numel(parcel_names));
for parcel_n = 1:numel(parcel_patterns)
    for sub_n = 1:size(cv_es, 1)
        subject         = cv_es{sub_n, 3};
        result_homo_dir = fullfile(results_dir, 'homo', 'es_pp', subject);
        homo_files      = dir(fullfile(result_homo_dir, ...
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
        sub_homo_es(sub_n, parcel_n) = mean(fold_homos, 'omitnan');
    end
end

p_vals_es = [];
pair_names_es = {};
for i = 1:numel(parcel_names)
    for j = (i+1):numel(parcel_names)
        [~, p]               = ttest(sub_homo_es(:, i), sub_homo_es(:, j));
        p_vals_es(end+1)     = p; %#ok<SAGROW>
        pair_names_es{end+1} = [parcel_names{i} ' vs ' parcel_names{j}]; %#ok<SAGROW>
    end
end
p_fdr_es = mafdr(p_vals_es, 'BHFDR', true);

fprintf('\n--- es_pp Mean FC Homogeneity per Condition ---\n');
for parcel_n = 1:numel(parcel_names)
    fprintf('%s: %.4f +/- %.4f\n', parcel_names{parcel_n}, ...
        nanmean(sub_homo_es(:, parcel_n)), nanstd(sub_homo_es(:, parcel_n)));
end

fprintf('\n--- es_pp Paired T-Test Results (FDR corrected, q<0.05) ---\n');
for k = 1:numel(pair_names_es)
    sig = '';
    if p_fdr_es(k) < 0.05
        sig = ' *';
    end
    fprintf('%s: p_fdr=%.4f%s\n', pair_names_es{k}, p_fdr_es(k), sig);
end

% Save es_pp results to CSV
summary_file = fullfile(results_dir, 'homo_es_pp_summary.csv');
fid = fopen(summary_file, 'w');
fprintf(fid, 'condition,mean_homo,sd_homo\n');
for p = 1:numel(parcel_names)
    fprintf(fid, '%s,%.6f,%.6f\n', parcel_names{p}, ...
        nanmean(sub_homo_es(:, p)), nanstd(sub_homo_es(:, p)));
end
fclose(fid);
fprintf('\nSaved: %s\n', summary_file);

stats_file = fullfile(results_dir, 'homo_es_pp_stats.csv');
fid = fopen(stats_file, 'w');
fprintf(fid, 'pair,p_uncorr,p_fdr,sig_uncorr,sig_fdr\n');
for k = 1:numel(pair_names_es)
    fprintf(fid, '%s,%.6f,%.6f,%d,%d\n', pair_names_es{k}, ...
        p_vals_es(k), p_fdr_es(k), p_vals_es(k) < 0.05, p_fdr_es(k) < 0.05);
end
fclose(fid);
fprintf('Saved: %s\n', stats_file);

% ====================================
%% Cleanup Paths
% ====================================

rmpath(proj_dir);
rmpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Local Functions
% ====================================

