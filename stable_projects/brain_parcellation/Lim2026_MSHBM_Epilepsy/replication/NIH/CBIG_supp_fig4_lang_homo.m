% ====================================
%% Help
% ====================================
%
% CBIG_supp_fig4_lang_homo
%
% This script replicates Supplementary Figure 4: FC homogeneity and
% inhomogeneity restricted to the Language A and Language B networks.
% Four panels are produced:
%   A - NIH test set resting-state homogeneity
%   B - esfMRI pre-operative resting-state homogeneity
%   C - esfMRI post-operative resting-state homogeneity
%   D - esfMRI post-operative electrical stimulation inhomogeneity
%
% For each panel, six conditions are compared (HCP group, HCP individual,
% Du group, Du individual, NIH group, NIH individual). Within-language-
% network homogeneity for Language A and B is computed separately then
% combined as a network-size-weighted average. Statistics use paired
% t-tests with BH-FDR correction (q<0.05). Each panel is displayed as
% a figure with a boxplot (left) and FDR-corrected p-value table (right).
%
% PREREQUISITES:
%   Run CBIG_test_MSHBM_Epilepsy.m to generate NIH parcellations and
%   homolist files.
%   Run CBIG_test_MSHBM_esfmri_homo.m to generate esfMRI parcellations
%   and homolist files.
%   Run CBIG_test_MSHBM_esfmri_glm.m and CBIG_test_MSHBM_esfmri_inhomo.m
%   to generate GLM gamma maps (for Panel D).
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   NIH parcellations and homolist files in:
%   $proj_dir/replication/NIH/input/{subject}/mshbm/{prior}_cv_{n}/
%   $proj_dir/replication/NIH/input/{subject}/homo/
%   esfMRI parcellations and homolist files in:
%   $proj_dir/replication/esfmri/input/{rs_pp|es_pp}/{subject}/mshbm/
%   $proj_dir/replication/esfmri/input/{rs_pp|es_pp}/{subject}/homo/
%   GLM gamma maps in:
%   $proj_dir/replication/esfmri/input/es_pp/{subject}/glm/
%   Group parcellation labels in $proj_dir/replication/mshbm/{HCP|Du|NIH}/
%
% OUTPUTS:
%   Prints mean +/- SD language homogeneity per condition and BH-FDR
%   corrected paired t-test results to console.
%   Displays Supplementary Figure 4 panels A-D on screen (not saved).
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

nih_input_dir  = fullfile(data_dir, 'NIH/input');
esfmri_dir     = fullfile(data_dir, 'esfmri/input');
mshbm_ref_dir  = fullfile(proj_dir, 'replication/mshbm');

addpath(fullfile(proj_dir, 'replication/NIH'));

% Group parcellation labels
grppar_paths = {fullfile(mshbm_ref_dir, 'HCP/group.mat'), ...
                fullfile(mshbm_ref_dir, 'Du/group_consensus_raw.mat'), ...
                fullfile(mshbm_ref_dir, 'NIH/group.mat')};

% Medial wall mask
load(fullfile(mshbm_ref_dir, 'union_medial_wall.mat'), 'mw_lh', 'mw_rh');

% Language network labels per prior (langA, langB): HCP, Du, NIH
lang_labels = {[5, 4]; [9, 3]; [8, 11]};

% mshbm_flags: folder name prefix for parcellations (esfmri uses HCP/Du/NIH;
% NIH test uses HCP_15net/Du_15net/NIH_34 due to how those were generated)
mshbm_flags     = {'HCP', 'Du', 'NIH'};
nih_parcel_flags = {'HCP_15net', 'Du_15net', 'NIH_34'};
parcel_names = {'HCP_group', 'HCP_IP', 'Du_group', 'Du_IP', 'NIH_group', 'NIH_IP'};
run_pattern  = 'bld([\d]+)';

% ====================================
%% Panel A: NIH Test Set Resting-State Homogeneity
% ====================================

fprintf('\n====== Panel A: NIH test set language homogeneity ======\n');

fid = fopen(fullfile(nih_input_dir, 'NIH_test_14_subid.csv'));
fgetl(fid);
subid_nih = textscan(fid, '%s', 'Delimiter', '\t');
fclose(fid);
subid_nih = subid_nih{1};

% Build leave-one-run-out CV structure
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
fprintf('CV subjects: %d\n', size(cv_nih, 1));

sub_homo_nih = compute_lang_homo(cv_nih, nih_input_dir, nih_input_dir, ...
    nih_parcel_flags, grppar_paths, lang_labels, mw_lh, mw_rh, run_pattern);

[p_mat_nih, q_mat_nih, dir_mat_nih] = compute_stats(sub_homo_nih);
print_results('Panel A (NIH test set)', sub_homo_nih, parcel_names, p_mat_nih, q_mat_nih);

% ====================================
%% Panel B: esfMRI Pre-op Resting-State Homogeneity
% ====================================

fprintf('\n====== Panel B: esfMRI pre-op language homogeneity ======\n');

rspp_dir = fullfile(esfmri_dir, 'rs_pp');
fid = fopen(fullfile(esfmri_dir, 'rsfmri_subid.csv'));
fgetl(fid);
subid_rs = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid_rs = subid_rs{1};

mc_rs = load(fullfile(rspp_dir, 'mc_test.mat'), 'mc_test');
mc_test_rs = mc_rs.mc_test;

subj_perm_rs = cell(numel(subid_rs), 3);
for sub_n = 1:numel(subid_rs)
    subject = subid_rs{sub_n};
    subj    = strrep(subject, '-', '_');
    if ~isfield(mc_test_rs, subj) || sum(mc_test_rs.(subj).pass) < 2
        continue
    end
    surf_files = dir(fullfile(rspp_dir, subject, 'surf', 'lh.*fs6_sm6_censored.nii.gz'));
    surf_names = {};
    for i = 1:numel(surf_files)
        if mc_test_rs.(subj).pass(i) == 0; continue; end
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
fprintf('CV subjects: %d\n', size(cv_rs, 1));

sub_homo_rs = compute_lang_homo(cv_rs, rspp_dir, rspp_dir, ...
    mshbm_flags, grppar_paths, lang_labels, mw_lh, mw_rh, run_pattern);

[p_mat_rs, q_mat_rs, dir_mat_rs] = compute_stats(sub_homo_rs);
print_results('Panel B (esfMRI pre-op)', sub_homo_rs, parcel_names, p_mat_rs, q_mat_rs);

% ====================================
%% Panel C: esfMRI Post-op Resting-State Homogeneity
% ====================================

fprintf('\n====== Panel C: esfMRI post-op language homogeneity ======\n');

espp_dir = fullfile(esfmri_dir, 'es_pp');
fid = fopen(fullfile(esfmri_dir, 'esfmri_subid.csv'));
fgetl(fid);
subid_es = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid_es = subid_es{1};

mc_es = load(fullfile(espp_dir, 'mc_test.mat'), 'mc_test');
mc_test_es = mc_es.mc_test;

subj_perm_es = cell(numel(subid_es), 3);
for sub_n = 1:numel(subid_es)
    subject = subid_es{sub_n};
    subj    = strrep(subject, '-', '_');
    if ~isfield(mc_test_es, subj) || sum(mc_test_es.(subj).pass) < 2
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
fprintf('CV subjects: %d\n', size(cv_es, 1));

sub_homo_es = compute_lang_homo(cv_es, espp_dir, espp_dir, ...
    mshbm_flags, grppar_paths, lang_labels, mw_lh, mw_rh, run_pattern);

[p_mat_es, q_mat_es, dir_mat_es] = compute_stats(sub_homo_es);
print_results('Panel C (esfMRI post-op)', sub_homo_es, parcel_names, p_mat_es, q_mat_es);

% ====================================
%% Panel D: esfMRI Post-op Electrical Stimulation Inhomogeneity
% ====================================

fprintf('\n====== Panel D: esfMRI post-op language inhomogeneity ======\n');

sub_inhomo_es = compute_lang_inhomo(cv_es, espp_dir, ...
    mshbm_flags, grppar_paths, lang_labels, mw_lh, mw_rh);

[p_mat_inhomo, q_mat_inhomo, dir_mat_inhomo] = compute_stats(sub_inhomo_es);
print_results('Panel D (esfMRI inhomo)', sub_inhomo_es, parcel_names, ...
    p_mat_inhomo, q_mat_inhomo);

% ====================================
%% Plot Supplementary Figure 4
% ====================================

panel_data    = {sub_homo_nih, sub_homo_rs, sub_homo_es, sub_inhomo_es};
panel_titles  = {'A: NIH test set rs-fMRI', 'B: esfMRI pre-op rs-fMRI', ...
                 'C: esfMRI post-op rs-fMRI', 'D: esfMRI post-op ES inhomogeneity'};
panel_ylabels = {'Language A+B Homogeneity', 'Language A+B Homogeneity', ...
                 'Language A+B Homogeneity', 'Language A+B Inhomogeneity'};
panel_higher  = {true, true, true, false};
panel_pmats   = {p_mat_nih,    p_mat_rs,    p_mat_es,    p_mat_inhomo};
panel_qmats   = {q_mat_nih,    q_mat_rs,    q_mat_es,    q_mat_inhomo};
panel_dirmats = {dir_mat_nih,  dir_mat_rs,  dir_mat_es,  dir_mat_inhomo};

cond_labels = {'HCP avg', 'HCP ind', 'Du avg', 'Du ind', 'NIH avg', 'NIH ind'};
model_labels = {'HCP', 'Du', 'NIH'};

for panel_n = 1:4
    data       = panel_data{panel_n};
    p_mat      = panel_pmats{panel_n};
    q_mat      = panel_qmats{panel_n};
    dir_mat    = panel_dirmats{panel_n};
    higher     = panel_higher{panel_n};
    n_subs     = size(data, 1);
    n_cond     = size(data, 2);

    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % --- Left: boxplot ---
    ax1 = subplot(1, 2, 1);
    hold(ax1, 'on');

    y_min = min(data(:), [], 'omitnan');
    y_max = max(data(:), [], 'omitnan');
    y_pad = 0.15 * (y_max - y_min);

    for m = 1:3
        grp_col = (m-1)*2 + 1;
        ind_col = (m-1)*2 + 2;
        for sub_n = 1:n_subs
            plot(ax1, [grp_col ind_col], [data(sub_n, grp_col) data(sub_n, ind_col)], ...
                '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        end
    end

    bw = 0.3;
    for c = 1:n_cond
        cdata = data(~isnan(data(:, c)), c);
        if isempty(cdata); continue; end
        q1 = quantile(cdata, 0.25); q2 = quantile(cdata, 0.50); q3 = quantile(cdata, 0.75);
        iqr_val = q3 - q1;
        w_lo = max(min(cdata), q1 - 1.5*iqr_val);
        w_hi = min(max(cdata), q3 + 1.5*iqr_val);
        rectangle(ax1, 'Position', [c-bw, q1, 2*bw, q3-q1], ...
            'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 1);
        line(ax1, [c-bw, c+bw], [q2, q2], 'Color', 'k', 'LineWidth', 1.5);
        line(ax1, [c, c], [w_lo, q1], 'Color', 'k', 'LineWidth', 1);
        line(ax1, [c, c], [q3, w_hi], 'Color', 'k', 'LineWidth', 1);
        line(ax1, [c-bw/2, c+bw/2], [w_lo, w_lo], 'Color', 'k', 'LineWidth', 1);
        line(ax1, [c-bw/2, c+bw/2], [w_hi, w_hi], 'Color', 'k', 'LineWidth', 1);
        x_jit = c + (rand(n_subs, 1) - 0.5) * 0.08;
        scatter(ax1, x_jit, data(:, c), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end

    for m = 1:3
        grp_col = (m-1)*2 + 1; ind_col = (m-1)*2 + 2;
        x_mid = (grp_col + ind_col) / 2;
        text(ax1, x_mid, y_min - 0.07*(y_max-y_min), model_labels{m}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14, ...
            'FontWeight', 'bold', 'FontName', 'Arial');
        line(ax1, [grp_col-0.35 ind_col+0.35], ...
            [y_min-0.04*(y_max-y_min), y_min-0.04*(y_max-y_min)], ...
            'Color', 'k', 'LineWidth', 1);
    end

    ax1.YLim = [y_min - y_pad, y_max + y_pad];
    ax1.XLim = [0.5, n_cond + 0.5];
    ax1.XTick = 1:n_cond;
    ax1.XTickLabel = {'avg', 'ind', 'avg', 'ind', 'avg', 'ind'};
    ax1.XAxis.FontSize = 14; ax1.YAxis.FontSize = 14;
    ax1.FontName = 'Arial'; ax1.Box = 'off'; ax1.TickDir = 'out';
    ylabel(ax1, panel_ylabels{panel_n}, 'FontName', 'Arial', 'FontSize', 14);
    title(ax1, panel_titles{panel_n}, 'FontName', 'Arial', 'FontSize', 14, ...
        'FontWeight', 'bold');

    % --- Right: p-value table ---
    ax2 = subplot(1, 2, 2);
    axis(ax2, 'off');
    hold(ax2, 'on');

    cell_w = 0.80 / n_cond;
    cell_h = 0.70 / n_cond;
    x_off  = 0.18;
    y_off  = 0.05;

    for j = 1:n_cond
        x_ctr = x_off + (j - 0.5) * cell_w;
        text(ax2, x_ctr, y_off + n_cond*cell_h + 0.04, cond_labels{j}, ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontName', 'Arial');
    end
    for m = 1:3
        grp_col = (m-1)*2 + 1; ind_col = (m-1)*2 + 2;
        x_mid = x_off + ((grp_col + ind_col)/2 - 0.5) * cell_w;
        text(ax2, x_mid, y_off + n_cond*cell_h + 0.10, model_labels{m}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14, ...
            'FontWeight', 'bold', 'FontName', 'Arial');
        line(ax2, [x_off+(grp_col-1)*cell_w+0.005, x_off+ind_col*cell_w-0.005], ...
            repmat(y_off + n_cond*cell_h + 0.07, 1, 2), 'Color', 'k', 'LineWidth', 1);
    end

    if higher; better_dir = 1; else; better_dir = -1; end

    for i = 1:n_cond
        y_ctr = y_off + (n_cond - i + 0.5) * cell_h;
        text(ax2, x_off - 0.01, y_ctr, cond_labels{i}, ...
            'HorizontalAlignment', 'right', 'FontSize', 12, 'FontName', 'Arial');
        for j = 1:n_cond
            x_pos = x_off + (j-1)*cell_w + 0.003;
            y_pos = y_off + (n_cond - i)*cell_h + 0.003;
            cw = cell_w - 0.006; ch = cell_h - 0.006;
            if i == j
                rectangle(ax2, 'Position', [x_pos y_pos cw ch], ...
                    'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'w', 'LineWidth', 1);
            else
                q_val = q_mat(i, j); dir_val = dir_mat(i, j);
                p_raw = p_mat(i, j);
                if isnan(p_raw)
                    p_str = '';
                elseif p_raw < 0.001
                    p_str = 'p<0.001';
                else
                    p_str = sprintf('p=%.3f', p_raw);
                end
                if isnan(q_val) || q_val >= 0.05
                    bg = [0.75 0.75 0.75];
                elseif dir_val * better_dir > 0
                    bg = [0.35 0.75 0.35];
                else
                    bg = [0.90 0.25 0.25];
                end
                rectangle(ax2, 'Position', [x_pos y_pos cw ch], ...
                    'FaceColor', bg, 'EdgeColor', 'w', 'LineWidth', 1);
                if ~isempty(p_str)
                    text(ax2, x_pos + cw/2, y_pos + ch/2, p_str, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                        'FontSize', 11, 'FontName', 'Arial');
                end
            end
        end
    end
    for m = 1:3
        grp_row = (m-1)*2 + 1; ind_row = (m-1)*2 + 2;
        y_mid = y_off + (n_cond - (grp_row + ind_row)/2) * cell_h;
        text(ax2, x_off - 0.12, y_mid, model_labels{m}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14, ...
            'FontWeight', 'bold', 'FontName', 'Arial', 'Rotation', 90);
    end
end

% ====================================
%% Cleanup
% ====================================

rmpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Local Functions
% ====================================

function sub_homo = compute_lang_homo(cv, surf_base, parcel_base, mshbm_flags, ...
    grppar_paths, lang_labels, mw_lh, mw_rh, run_pattern)
% Compute language A+B weighted homogeneity per subject per condition.

n_subs   = size(cv, 1);
sub_homo = NaN(n_subs, 6);

for sub_n = 1:n_subs
    subject   = cv{sub_n, 3};
    num_folds = numel(cv{sub_n, 2});
    homo_dir  = fullfile(surf_base, subject, 'homo');

    fold_homo = NaN(num_folds, 6);

    % Pre-load group parcellations
    grp_labels = cell(numel(mshbm_flags), 2);
    for pg = 1:numel(mshbm_flags)
        gp = load(grppar_paths{pg});
        lh_lbl = double(gp.lh_labels); lh_lbl(mw_lh) = 0;
        rh_lbl = double(gp.rh_labels); rh_lbl(mw_rh) = 0;
        grp_labels{pg, 1} = lh_lbl;
        grp_labels{pg, 2} = rh_lbl;
    end

    for fold_n = 1:num_folds
        test_file = cv{sub_n, 2}{fold_n}{1};
        bld       = regexp(test_file, run_pattern, 'match');
        lhlist    = fullfile(homo_dir, ['lh_homolist_' bld{1} '.txt']);
        rhlist    = fullfile(homo_dir, ['rh_homolist_' bld{1} '.txt']);

        if ~exist(lhlist, 'file'); continue; end

        for pm = 1:numel(mshbm_flags)
            langA = lang_labels{pm}(1); langB = lang_labels{pm}(2);
            grp_col = (pm-1)*2 + 1; ind_col = pm*2;

            % Group parcellation homogeneity
            fold_homo(fold_n, grp_col) = lang_homo_AB( ...
                grp_labels{pm, 1}, grp_labels{pm, 2}, lhlist, rhlist, langA, langB);

            % Individual parcellation homogeneity
            parcel_file = fullfile(parcel_base, subject, 'mshbm', ...
                [mshbm_flags{pm} '_cv_' num2str(fold_n)], ...
                'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');
            if ~exist(parcel_file, 'file'); continue; end
            ip = load(parcel_file);
            lh_lbl = double(ip.lh_labels); lh_lbl(mw_lh) = 0;
            rh_lbl = double(ip.rh_labels); rh_lbl(mw_rh) = 0;
            fold_homo(fold_n, ind_col) = lang_homo_AB( ...
                lh_lbl, rh_lbl, lhlist, rhlist, langA, langB);
        end

        fprintf('  %s fold %d done\n', subject, fold_n);
    end

    sub_homo(sub_n, :) = mean(fold_homo, 1, 'omitnan');
end
end

function sub_inhomo = compute_lang_inhomo(cv_es, espp_dir, mshbm_flags, ...
    grppar_paths, lang_labels, mw_lh, mw_rh)
% Compute language A+B weighted inhomogeneity per subject per condition.

n_subs     = size(cv_es, 1);
sub_inhomo = NaN(n_subs, 6);

for sub_n = 1:n_subs
    subject   = cv_es{sub_n, 3};
    num_folds = numel(cv_es{sub_n, 2});
    glm_dir   = fullfile(espp_dir, subject, 'glm');

    glm_lh = dir(fullfile(glm_dir, '*_glmfit_flob3_lh'));
    glm_rh = dir(fullfile(glm_dir, '*_glmfit_flob3_rh'));

    if numel(glm_lh) < num_folds
        fprintf('  Fewer GLM outputs than folds for %s; skipping.\n', subject);
        continue
    end

    grp_labels = cell(numel(mshbm_flags), 2);
    for pg = 1:numel(mshbm_flags)
        gp = load(grppar_paths{pg});
        lh_lbl = double(gp.lh_labels); lh_lbl(mw_lh) = 0;
        rh_lbl = double(gp.rh_labels); rh_lbl(mw_rh) = 0;
        grp_labels{pg, 1} = lh_lbl;
        grp_labels{pg, 2} = rh_lbl;
    end

    fold_inhomo = NaN(num_folds, 6);

    for fold_n = 1:num_folds
        lh_gamma = fullfile(glm_lh(fold_n).folder, glm_lh(fold_n).name, ...
            'stim_contrast_pw', 'gamma.nii.gz');
        rh_gamma = fullfile(glm_rh(fold_n).folder, glm_rh(fold_n).name, ...
            'stim_contrast_pw', 'gamma.nii.gz');
        if ~exist(lh_gamma, 'file') || ~exist(rh_gamma, 'file'); continue; end

        lh_s = MRIread(lh_gamma); rh_s = MRIread(rh_gamma);
        lh_GLM = double(lh_s.vol(:)); rh_GLM = double(rh_s.vol(:));

        for pm = 1:numel(mshbm_flags)
            langA = lang_labels{pm}(1); langB = lang_labels{pm}(2);
            grp_col = (pm-1)*2 + 1; ind_col = pm*2;

            fold_inhomo(fold_n, grp_col) = lang_inhomo_AB( ...
                grp_labels{pm, 1}, grp_labels{pm, 2}, lh_GLM, rh_GLM, langA, langB);

            parcel_file = fullfile(espp_dir, subject, 'mshbm', ...
                [mshbm_flags{pm} '_cv_' num2str(fold_n)], ...
                'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');
            if ~exist(parcel_file, 'file'); continue; end
            ip = load(parcel_file);
            lh_lbl = double(ip.lh_labels); lh_lbl(mw_lh) = 0;
            rh_lbl = double(ip.rh_labels); rh_lbl(mw_rh) = 0;
            fold_inhomo(fold_n, ind_col) = lang_inhomo_AB( ...
                lh_lbl, rh_lbl, lh_GLM, rh_GLM, langA, langB);
        end

        fprintf('  %s fold %d done\n', subject, fold_n);
    end

    sub_inhomo(sub_n, :) = mean(fold_inhomo, 1, 'omitnan');
end
end

function homo_AB = lang_homo_AB(lh_labels, rh_labels, lhlist, rhlist, langA, langB)
% Weighted language A+B homogeneity for one fold.
lh_A = make_lang_mask(lh_labels, langA); rh_A = make_lang_mask(rh_labels, langA);
lh_B = make_lang_mask(lh_labels, langB); rh_B = make_lang_mask(rh_labels, langB);
n_A = sum(lh_A(:) > 0) + sum(rh_A(:) > 0);
n_B = sum(lh_B(:) > 0) + sum(rh_B(:) > 0);
if n_A + n_B == 0; homo_AB = NaN; return; end
h_A = CBIG_ParcellationHomogeneity_FS_meantimecourse(lh_A, rh_A, 'fsaverage6', lhlist, rhlist);
h_B = CBIG_ParcellationHomogeneity_FS_meantimecourse(lh_B, rh_B, 'fsaverage6', lhlist, rhlist);
homo_AB = (n_A * h_A + n_B * h_B) / (n_A + n_B);
end

function inhomo_AB = lang_inhomo_AB(lh_labels, rh_labels, lh_GLM, rh_GLM, langA, langB)
% Weighted language A+B inhomogeneity for one fold.
[inh_A, n_A] = lang_net_inhomo(lh_labels, rh_labels, langA, lh_GLM, rh_GLM);
[inh_B, n_B] = lang_net_inhomo(lh_labels, rh_labels, langB, lh_GLM, rh_GLM);
if n_A + n_B == 0; inhomo_AB = NaN; return; end
inhomo_AB = (n_A * inh_A + n_B * inh_B) / (n_A + n_B);
end

function [inhomo, n] = lang_net_inhomo(lh_labels, rh_labels, label, lh_GLM, rh_GLM)
% Within-network std of z-scored GLM map, weighted by network size.
wb_mean = mean([lh_GLM(:); rh_GLM(:)], 'omitnan');
wb_std  = std([lh_GLM(:); rh_GLM(:)], 'omitnan');
lh_z = (lh_GLM - wb_mean) / wb_std; rh_z = (rh_GLM - wb_mean) / wb_std;
lh_mask = (lh_labels == label); rh_mask = (rh_labels == label);
n = sum(lh_mask) + sum(rh_mask);
if n == 0; inhomo = NaN; return; end
all_vals = [lh_z(lh_mask); rh_z(rh_mask)];
net_std  = std(all_vals, 'omitnan');
net_wt   = n / (sum(lh_labels > 0) + sum(rh_labels > 0));
inhomo   = net_std * net_wt;
end

function mask = make_lang_mask(labels, label)
% Return binary mask with specified label set to 1, all others to 0.
mask = double(labels == label);
end

function [p_mat, q_mat, dir_mat] = compute_stats(data)
% Pairwise t-tests and BH-FDR correction.
n_cond  = size(data, 2);
p_mat   = NaN(n_cond, n_cond);
dir_mat = zeros(n_cond, n_cond);
for i = 1:n_cond
    for j = 1:n_cond
        if i == j; continue; end
        [~, p, ~, stats] = ttest(data(:, i), data(:, j));
        p_mat(i, j)   = p;
        dir_mat(i, j) = sign(stats.tstat);
    end
end
mask      = ~isnan(p_mat) & ~eye(n_cond, 'logical');
p_vec     = p_mat(mask);
q_vec     = mafdr(p_vec, 'BHFDR', true);
q_mat     = NaN(n_cond, n_cond);
q_mat(mask) = q_vec;
end

function print_results(label, data, parcel_names, p_mat, q_mat)
fprintf('\n--- %s: Mean +/- SD ---\n', label);
for p = 1:numel(parcel_names)
    fprintf('%s: %.4f +/- %.4f\n', parcel_names{p}, ...
        mean(data(:, p), 'omitnan'), std(data(:, p), 'omitnan'));
end
fprintf('--- FDR-corrected p-values (significant pairs) ---\n');
for i = 1:numel(parcel_names)
    for j = (i+1):numel(parcel_names)
        if ~isnan(q_mat(i, j)) && q_mat(i, j) < 0.05
            fprintf('  %s vs %s: p_fdr=%.4f\n', ...
                parcel_names{i}, parcel_names{j}, q_mat(i, j));
        end
    end
end
end
