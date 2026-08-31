% ====================================
%% Help
% ====================================
%
% CBIG_supp_fig5_model_compare
%
% This script replicates Supplementary Figure 5: FC homogeneity and
% inhomogeneity of the NIH individual-specific networks as a function of
% the number of training subjects used to train MS-HBM (5, 10, 15, 20, 25,
% 30, 34). Four panels are plotted in a single figure:
%   A - NIH test set resting-state homogeneity
%   B - esfMRI pre-operative resting-state homogeneity
%   C - esfMRI post-operative resting-state homogeneity
%   D - esfMRI post-operative electrical stimulation inhomogeneity
%
% Each panel shows mean +/- SEM across subjects as a line plot against
% training set size. Asterisks mark models that differ significantly from
% the n=34 reference after pooled BH-FDR correction (q<0.05) across all
% panels.
%
% The script first generates individual-specific parcellations for each
% model and computes FC homogeneity / inhomogeneity, then aggregates and
% plots the results.
%
% PREREQUISITES:
%   Review model priors must be present in:
%   replication/NIH/results/model_compare/{model}/
%     generate_individual_parcellations/priors/Params_Final.mat
%   Reference model (NIH_34) prior in:
%   replication/mshbm/NIH/Params_Final.mat
%   Preprocessed NIH fMRI surface files, motion censor files, and
%   homolist files in replication/NIH/input/{subject}/
%   esfMRI surface files, homolist files, and GLM gamma files in:
%   replication/esfmri/input/rs_pp/{subject}/ and es_pp/{subject}/
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%
% OUTPUTS:
%   Per-fold homo results (variable: homo_with_weight):
%   replication/NIH/results/model_compare/NIH_test_homo/{model}/{subject}/
%   replication/NIH/results/model_compare/preop_homo/{model}/{subject}/
%   replication/NIH/results/model_compare/postop_homo/{model}/{subject}/
%   Per-fold inhomo results (variable: network_std_mean):
%   replication/NIH/results/model_compare/postop_inhomo/{model}/{subject}/
%   Prints mean +/- SD per model per panel and FDR-corrected p-values to
%   console. Displays Supplementary Figure 5 on screen (not saved).
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Last tested on 20 May 2026.

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, ...
    'stable_projects/brain_parcellation/Lim2026_MSHBM_Epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');

nih_input_dir = fullfile(data_dir, 'NIH/input');
esfmri_dir    = fullfile(data_dir, 'esfmri/input');
mshbm_ref_dir = fullfile(proj_dir, 'replication/mshbm');
mc_dir        = fullfile(proj_dir, 'replication/NIH/results/model_compare');
rspp_dir      = fullfile(esfmri_dir, 'rs_pp');
espp_dir      = fullfile(esfmri_dir, 'es_pp');

medial_wall_path = fullfile(mshbm_ref_dir, 'union_medial_wall.mat');

addpath(proj_dir);
addpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Model Configurations
% ====================================

model_names = {'review_5', 'review_10', 'review_15', 'review_20', ...
               'review_25', 'review_30', 'NIH_34'};
model_sizes = [5, 10, 15, 20, 25, 30, 34];
n_models    = numel(model_names);
ref_idx     = n_models;    % n=34 is the reference

% Prior path for each model (review_X from mc_dir; reference from mshbm_ref_dir)
prior_paths = cell(1, n_models);
for m = 1:n_models - 1
    prior_paths{m} = fullfile(mc_dir, model_names{m}, ...
        'generate_individual_parcellations', 'priors', 'Params_Final.mat');
end
prior_paths{n_models} = fullfile(mshbm_ref_dir, 'NIH', 'Params_Final.mat');

w = '80';
c = '10';
run_pattern = 'bld([\d]+)';

% ====================================
%% Build NIH Test CV Structure
% ====================================

fid = fopen(fullfile(nih_input_dir, 'NIH_test_14_subid.csv'));
fgetl(fid);
subid_nih = textscan(fid, '%s', 'Delimiter', '\t');
fclose(fid);
subid_nih = subid_nih{1};

subj_perm_nih = cell(numel(subid_nih), 3);
for sub_n = 1:numel(subid_nih)
    subject = subid_nih{sub_n};
    surf_files = dir(fullfile(nih_input_dir, subject, 'surf', 'lh*sm6.nii.gz'));
    surf_names = {};
    for i = 1:numel(surf_files)
        surf_names{i} = surf_files(i).name(4:end);
    end
    subj_perm_nih{sub_n, 1} = cellfun(@(f) surf_names(~strcmp(surf_names, f)), ...
        surf_names, 'UniformOutput', false);
    subj_perm_nih{sub_n, 2} = cellfun(@(f) surf_names(strcmp(surf_names, f)), ...
        surf_names, 'UniformOutput', false);
    subj_perm_nih{sub_n, 3} = subject;
end
rows_rm = cellfun(@(x) numel(x) <= 1, subj_perm_nih(:, 1)) | ...
    cellfun(@(x) numel(x) <= 1, subj_perm_nih(:, 2));
cv_nih = subj_perm_nih(~rows_rm, :);

% ====================================
%% Build esfMRI Pre-op CV Structure
% ====================================

fid = fopen(fullfile(esfmri_dir, 'rsfmri_subid.csv'));
fgetl(fid);
subid_rs = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid_rs = subid_rs{1};
mc_rs = load(fullfile(rspp_dir, 'mc_test.mat'), 'mc_test');

subj_perm_rs = cell(numel(subid_rs), 3);
for sub_n = 1:numel(subid_rs)
    subject = subid_rs{sub_n};
    subj    = strrep(subject, '-', '_');
    if ~isfield(mc_rs.mc_test, subj) || sum(mc_rs.mc_test.(subj).pass) < 2
        continue
    end
    surf_files = dir(fullfile(rspp_dir, subject, 'surf', 'lh.*fs6_sm6_censored.nii.gz'));
    surf_names = {};
    for i = 1:numel(surf_files)
        if mc_rs.mc_test.(subj).pass(i) == 0; continue; end
        surf_names{end+1} = surf_files(i).name(4:end); %#ok<SAGROW>
    end
    subj_perm_rs{sub_n, 1} = cellfun(@(f) surf_names(~strcmp(surf_names, f)), ...
        surf_names, 'UniformOutput', false);
    subj_perm_rs{sub_n, 2} = cellfun(@(f) surf_names(strcmp(surf_names, f)), ...
        surf_names, 'UniformOutput', false);
    subj_perm_rs{sub_n, 3} = subject;
end
rows_rm = cellfun(@(x) numel(x) <= 1, subj_perm_rs(:, 1)) | ...
    cellfun(@(x) numel(x) <= 1, subj_perm_rs(:, 2));
cv_rs = subj_perm_rs(~rows_rm, :);

% ====================================
%% Build esfMRI Post-op CV Structure
% ====================================

fid = fopen(fullfile(esfmri_dir, 'esfmri_subid.csv'));
fgetl(fid);
subid_es = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid_es = subid_es{1};
mc_es = load(fullfile(espp_dir, 'mc_test.mat'), 'mc_test');

subj_perm_es = cell(numel(subid_es), 3);
for sub_n = 1:numel(subid_es)
    subject = subid_es{sub_n};
    subj    = strrep(subject, '-', '_');
    if ~isfield(mc_es.mc_test, subj) || sum(mc_es.mc_test.(subj).pass) < 2
        continue
    end
    surf_files = dir(fullfile(espp_dir, subject, 'surf', 'lh.*_nostim_discard.nii.gz'));
    surf_names = {};
    for i = 1:numel(surf_files)
        surf_names{end+1} = surf_files(i).name(4:end); %#ok<SAGROW>
    end
    subj_perm_es{sub_n, 1} = cellfun(@(f) surf_names(~strcmp(surf_names, f)), ...
        surf_names, 'UniformOutput', false);
    subj_perm_es{sub_n, 2} = cellfun(@(f) surf_names(strcmp(surf_names, f)), ...
        surf_names, 'UniformOutput', false);
    subj_perm_es{sub_n, 3} = subject;
end
rows_rm = cellfun(@(x) numel(x) <= 1, subj_perm_es(:, 1)) | ...
    cellfun(@(x) numel(x) <= 1, subj_perm_es(:, 2));
cv_es = subj_perm_es(~rows_rm, :);

% ====================================
%% Load Medial Wall Mask
% ====================================

load(medial_wall_path, 'mw_lh', 'mw_rh');

% ====================================
%% Compute NIH Test Homogeneity
% ====================================

fprintf('\n====== Computing NIH test homogeneity ======\n');
for m = 1:n_models
    out_dir = fullfile(mc_dir, 'NIH_test_homo', model_names{m});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    compute_homo_section(cv_nih, nih_input_dir, out_dir, model_names{m}, ...
        prior_paths{m}, w, c, mw_lh, mw_rh, run_pattern, proj_dir, true);
end

% ====================================
%% Compute Pre-op Homogeneity
% ====================================

fprintf('\n====== Computing pre-op homogeneity ======\n');
for m = 1:n_models
    out_dir = fullfile(mc_dir, 'preop_homo', model_names{m});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    compute_homo_section(cv_rs, rspp_dir, out_dir, model_names{m}, ...
        prior_paths{m}, w, c, mw_lh, mw_rh, run_pattern, proj_dir, false);
end

% ====================================
%% Compute Post-op Homogeneity
% ====================================

fprintf('\n====== Computing post-op homogeneity ======\n');
for m = 1:n_models
    out_dir = fullfile(mc_dir, 'postop_homo', model_names{m});
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    compute_homo_section(cv_es, espp_dir, out_dir, model_names{m}, ...
        prior_paths{m}, w, c, mw_lh, mw_rh, run_pattern, proj_dir, false);
end

% ====================================
%% Compute Post-op Inhomogeneity
% ====================================

fprintf('\n====== Computing post-op inhomogeneity ======\n');
for m = 1:n_models
    model      = model_names{m};
    out_subdir = fullfile(mc_dir, 'postop_inhomo', model);
    if ~exist(out_subdir, 'dir'); mkdir(out_subdir); end

    for sub_n = 1:size(cv_es, 1)
        subject = cv_es{sub_n, 3};
        n_folds = numel(cv_es{sub_n, 1});
        sub_out = fullfile(out_subdir, subject);
        if ~exist(sub_out, 'dir'); mkdir(sub_out); end
        glm_dir = fullfile(espp_dir, subject, 'glm');

        for fold_n = 1:n_folds
            out_file = fullfile(sub_out, ...
                [subject '_inhomo_' model '_cv_' num2str(fold_n) '_network_stat.mat']);
            if exist(out_file, 'file')
                fprintf('Skip (exists): %s %s inhomo cv%d\n', subject, model, fold_n);
                continue
            end

            % Parcellation must exist from post-op homo step above
            parcel_path = fullfile(espp_dir, subject, 'mshbm', ...
                [model '_cv_' num2str(fold_n)], ...
                'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');
            if ~exist(parcel_path, 'file')
                fprintf('Parcellation not found: %s %s cv%d\n', subject, model, fold_n);
                continue
            end

            sub_IP    = load(parcel_path);
            lh_labels = double(sub_IP.lh_labels); lh_labels(mw_lh) = 0;
            rh_labels = double(sub_IP.rh_labels); rh_labels(mw_rh) = 0;

            % Match GLM file to the held-out run using bld number
            test_file = cv_es{sub_n, 2}{fold_n}{1};
            bld       = regexp(test_file, run_pattern, 'match');
            bld_str   = bld{1};
            lh_glm = dir(fullfile(glm_dir, ...
                ['*' bld_str '*_flob3_lh/stim_contrast_pw/gamma.nii.gz']));
            rh_glm = dir(fullfile(glm_dir, ...
                ['*' bld_str '*_flob3_rh/stim_contrast_pw/gamma.nii.gz']));
            if isempty(lh_glm) || isempty(rh_glm)
                fprintf('GLM not found: %s %s cv%d\n', subject, model, fold_n);
                continue
            end

            lh_GLM_s = MRIread(fullfile(lh_glm(1).folder, lh_glm(1).name));
            rh_GLM_s = MRIread(fullfile(rh_glm(1).folder, rh_glm(1).name));
            lh_GLM   = reshape(lh_GLM_s.vol, numel(lh_GLM_s.vol), 1);
            rh_GLM   = reshape(rh_GLM_s.vol, numel(rh_GLM_s.vol), 1);

            network_std_mean = compute_inhomo(lh_labels, rh_labels, lh_GLM, rh_GLM);
            save(out_file, 'network_std_mean');
            fprintf('Inhomo done: %s %s cv%d\n', subject, model, fold_n);
        end
    end
end

% ====================================
%% Aggregate Per-Subject Means
% ====================================

fprintf('\n====== Aggregating model comparison results ======\n');

nih_data    = aggregate_homo(cv_nih, model_names, mc_dir, 'NIH_test_homo',  'homo');
preop_data  = aggregate_homo(cv_rs,  model_names, mc_dir, 'preop_homo',     'homo');
postop_data = aggregate_homo(cv_es,  model_names, mc_dir, 'postop_homo',    'homo');
inhomo_data = aggregate_homo(cv_es,  model_names, mc_dir, 'postop_inhomo',  'inhomo');

% ====================================
%% Pooled FDR Correction
% ====================================

% All comparisons vs reference (n=34) pooled across all 4 panels
all_p = NaN(n_models-1, 4);
for m = 1:n_models-1
    [~, all_p(m, 1)] = ttest(nih_data(:, m),    nih_data(:, ref_idx));
    [~, all_p(m, 2)] = ttest(preop_data(:, m),  preop_data(:, ref_idx));
    [~, all_p(m, 3)] = ttest(postop_data(:, m), postop_data(:, ref_idx));
    [~, all_p(m, 4)] = ttest(inhomo_data(:, m), inhomo_data(:, ref_idx));
end

q_all = mafdr(all_p(:), 'BHFDR', true);
q_mat = reshape(q_all, n_models-1, 4);   % rows = models, cols = panels

fprintf('\n--- FDR-corrected q-values (each model vs n=34 reference) ---\n');
panel_names = {'NIH_test', 'preop_homo', 'postop_homo', 'postop_inhomo'};
for m = 1:n_models-1
    for p = 1:4
        fprintf('  %s vs 34  [%s]: p=%.4f  q=%.4f%s\n', ...
            model_names{m}, panel_names{p}, all_p(m, p), q_mat(m, p), ...
            char(double(q_mat(m, p) < 0.05) * ('*' - 0) + 0));
    end
end

% ====================================
%% Plot Supplementary Figure 5
% ====================================

panel_data    = {nih_data, preop_data, postop_data, inhomo_data};
panel_titles  = {'A: NIH test set', 'B: esfMRI pre-op rs-fMRI', ...
                 'C: esfMRI post-op rs-fMRI', 'D: esfMRI post-op stimulation'};
panel_ylabels = {'Homogeneity', 'Homogeneity', 'Homogeneity', 'Inhomogeneity'};

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

for panel_n = 1:4
    data   = panel_data{panel_n};
    q_vals = [q_mat(:, panel_n); NaN];   % NaN for reference (no asterisk)
    n_subs = size(data, 1);

    ax = subplot(2, 2, panel_n);
    hold(ax, 'on');

    grp_mean = mean(data, 1, 'omitnan');
    grp_sem  = std(data, 0, 1, 'omitnan') / sqrt(n_subs);

    errorbar(ax, model_sizes, grp_mean, grp_sem, 'ko-', ...
        'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');

    y_range = max(data(:), [], 'omitnan') - min(data(:), [], 'omitnan');
    y_top   = max(grp_mean + grp_sem, [], 'omitnan') + 0.05 * y_range;
    y_max   = y_top + 0.10 * y_range;
    y_min   = min(data(:), [], 'omitnan') - 0.05 * y_range;
    ax.YLim = [y_min, y_max];
    ax.XLim = [3, 36];

    for m = 1:n_models-1
        if ~isnan(q_vals(m)) && q_vals(m) < 0.05
            text(ax, model_sizes(m), y_top, '*', ...
                'HorizontalAlignment', 'center', 'FontSize', 16, ...
                'FontName', 'Arial', 'FontWeight', 'bold');
        end
    end

    xticks(ax, model_sizes);
    xticklabels(ax, arrayfun(@num2str, model_sizes, 'UniformOutput', false));
    xlabel(ax, 'Number of training subjects', 'FontName', 'Arial', 'FontSize', 14);
    ylabel(ax, panel_ylabels{panel_n}, 'FontName', 'Arial', 'FontSize', 14);
    title(ax, panel_titles{panel_n}, 'FontName', 'Arial', 'FontSize', 14, ...
        'FontWeight', 'bold');

    ax.FontName = 'Arial'; ax.FontSize = 14;
    ax.Box      = 'off';
    ax.TickDir  = 'out';
    ax.XGrid    = 'off'; ax.YGrid = 'off';
end

% ====================================
%% Cleanup Paths
% ====================================

rmpath(proj_dir);
rmpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Local Functions
% ====================================

function compute_homo_section(cv, data_dir, out_subdir, model, prior_path, ...
    w, c, mw_lh, mw_rh, run_pattern, proj_dir, use_censor)
% Generates parcellations and computes homo for all subjects in cv.
for sub_n = 1:size(cv, 1)
    subject  = cv{sub_n, 3};
    n_folds  = numel(cv{sub_n, 1});
    sub_out  = fullfile(out_subdir, subject);
    if ~exist(sub_out, 'dir'); mkdir(sub_out); end
    homo_dir = fullfile(data_dir, subject, 'homo');

    for fold_n = 1:n_folds
        out_file = fullfile(sub_out, ...
            [subject '_homo_' model '_cv_' num2str(fold_n) '.mat']);
        if exist(out_file, 'file')
            fprintf('Skip (exists): %s %s cv%d\n', subject, model, fold_n);
            continue
        end

        parcel_path = fullfile(data_dir, subject, 'mshbm', ...
            [model '_cv_' num2str(fold_n)], ...
            'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

        if ~exist(parcel_path, 'file')
            train_files = cv{sub_n, 1}{fold_n};
            lh_list = cellfun(@(f) fullfile(data_dir, subject, 'surf', ['lh.' f]), ...
                train_files, 'UniformOutput', false);
            rh_list = cellfun(@(f) fullfile(data_dir, subject, 'surf', ['rh.' f]), ...
                train_files, 'UniformOutput', false);
            if use_censor
                censor_list = cellfun(@(f) ...
                    get_censor_path(data_dir, subject, f, run_pattern), ...
                    train_files, 'UniformOutput', false);
            else
                censor_list = 'NONE';
            end
            out_dir = fullfile(data_dir, subject, 'mshbm', ...
                [model '_cv_' num2str(fold_n)]);
            run_parcellation(out_dir, prior_path, lh_list, rh_list, censor_list, ...
                w, c, proj_dir);
        end

        if ~exist(parcel_path, 'file')
            fprintf('Parcellation not found: %s %s cv%d\n', subject, model, fold_n);
            continue
        end

        sub_IP    = load(parcel_path);
        lh_labels = double(sub_IP.lh_labels); lh_labels(mw_lh) = 0;
        rh_labels = double(sub_IP.rh_labels); rh_labels(mw_rh) = 0;

        test_file = cv{sub_n, 2}{fold_n}{1};
        bld       = regexp(test_file, run_pattern, 'match');
        lhlist    = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
        rhlist    = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);
        if ~exist(lhlist, 'file')
            fprintf('Homolist not found: %s\n', lhlist); continue
        end

        homo_with_weight = CBIG_ParcellationHomogeneity_FS_meantimecourse(...
            lh_labels, rh_labels, 'fsaverage6', lhlist, rhlist);
        save(out_file, 'homo_with_weight');
        fprintf('Homo done: %s %s cv%d\n', subject, model, fold_n);
    end
end
end

function run_parcellation(out_dir, prior_path, lh_list, rh_list, censor_list, ...
    w, c, proj_dir)
% Calls CBIG_MSHBM_Epilepsy_LI with given parameters.
params             = struct();
params.project_dir = out_dir;
params.group_prior = prior_path;
params.lh_fMRI_list = lh_list;
params.rh_fMRI_list = rh_list;
params.censor_list  = censor_list;
params.target_mesh  = 'fsaverage6';
params.w            = w;
params.c            = c;
CBIG_MSHBM_Epilepsy_LI(params);
end

function network_std_mean = compute_inhomo(lh_labels, rh_labels, lh_GLM, rh_GLM)
% Computes weighted inhomogeneity: per-network std of z-scored GLM, weighted by network size.
wb_mean = mean([lh_GLM(:); rh_GLM(:)], 'omitnan');
wb_std  = std([lh_GLM(:);  rh_GLM(:)], 'omitnan');
if wb_std == 0
    network_std_mean = NaN; return
end
lh_z = (lh_GLM - wb_mean) / wb_std;
rh_z = (rh_GLM - wb_mean) / wb_std;

lh_net_std = zeros(1, 15);
rh_net_std = zeros(1, 15);
lh_net_sz  = zeros(1, 15);
rh_net_sz  = zeros(1, 15);
for net = 1:15
    lh_mask = lh_labels == net;
    rh_mask = rh_labels == net;
    lh_net_std(net) = std(lh_z(lh_mask), 'omitnan');
    rh_net_std(net) = std(rh_z(rh_mask), 'omitnan');
    lh_net_sz(net)  = sum(lh_mask);
    rh_net_sz(net)  = sum(rh_mask);
end
lh_wt = lh_net_sz / sum(lh_net_sz);
rh_wt = rh_net_sz / sum(rh_net_sz);
network_std_mean = mean(abs([lh_net_std .* lh_wt, rh_net_std .* rh_wt]));
end

function data = aggregate_homo(cv, model_names, mc_dir, subdir, mode)
% Loads per-fold homo or inhomo files, averages per subject, returns n_subs x n_models.
n_subs   = size(cv, 1);
n_models = numel(model_names);
data     = NaN(n_subs, n_models);

for m = 1:n_models
    model = model_names{m};
    for sub_n = 1:n_subs
        subject   = cv{sub_n, 3};
        n_folds   = numel(cv{sub_n, 1});
        fold_vals = NaN(n_folds, 1);

        for fold_n = 1:n_folds
            if strcmp(mode, 'inhomo')
                f = fullfile(mc_dir, subdir, model, subject, ...
                    [subject '_inhomo_' model '_cv_' num2str(fold_n) '_network_stat.mat']);
                if exist(f, 'file')
                    d = load(f, 'network_std_mean');
                    fold_vals(fold_n) = d.network_std_mean;
                end
            else
                f = fullfile(mc_dir, subdir, model, subject, ...
                    [subject '_homo_' model '_cv_' num2str(fold_n) '.mat']);
                if exist(f, 'file')
                    d = load(f, 'homo_with_weight');
                    fold_vals(fold_n) = d.homo_with_weight;
                end
            end
        end

        data(sub_n, m) = mean(fold_vals, 'omitnan');
    end
    fprintf('  Loaded %s / %s\n', subdir, model);
end
end

function censor_path = get_censor_path(data_dir, subject, surf_file, run_pattern)
% Returns motion censor file path for a given NIH surface fMRI filename.
bld = regexp(surf_file, run_pattern, 'match');
censor_path = fullfile(data_dir, subject, 'qc', ...
    [subject '_' bld{1} '_FDRMS0.2_DVARS50_motion_outliers.txt']);
end
