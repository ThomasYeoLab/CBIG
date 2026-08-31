% ====================================
%% Help
% ====================================
%
% CBIG_language_prediction
%
% This script generates full (all-run) individual-specific brain parcellations
% for NIH subjects with known language task laterality (L/B/R), then computes
% the language laterality index (LI) for each subject and parcellation scheme.
% LI is the hemispheric asymmetry in the number of vertices assigned to the
% combined language networks (labels A + B). One-way ANOVA with BH-FDR
% correction tests whether LI differs across task groups. Boxplots of LI
% and ROC curves for language laterality classification replicate Figure 5
% of Lim 2026.
% Last tested on 20 May 2026.
%
% PREREQUISITES:
%   Run CBIG_train_MSHBM_Epilepsy.m first to generate the NIH group prior in
%   replication/mshbm/NIH/.
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   Language laterality template (subid, task, n_runs):
%   $proj_dir/replication/NIH/input/NIH_lang_43_subid.csv
%   Preprocessed surface fMRI files (*sm6.nii.gz) in:
%   $proj_dir/replication/NIH/input/{subject}/surf/
%   Motion censor files in:
%   $proj_dir/replication/NIH/input/{subject}/qc/
%   Group priors (Params_Final.mat) in $proj_dir/replication/mshbm/{HCP|Du|NIH}/.
%
% OUTPUTS:
%   Individual parcellations (full, all runs, patient data — not released):
%   $proj_dir/replication/NIH/input/{subject}/mshbm/{HCP|Du|NIH}/
%   LI results table (subid, task, n_runs, {HCP|Du|NIH}_lang):
%   $proj_dir/replication/NIH/results/lang/NIH_language_prediction_results.xlsx
%   Prints ANOVA statistics (uncorrected and BH-FDR corrected) to console.
%   Displays boxplots of LI by task group and ROC curves (not saved to disk).
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, ...
    'stable_projects/brain_parcellation/Lim2026_MSHBM_Epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');
analysis_dir  = fullfile(data_dir, 'NIH/input');
results_dir   = fullfile(proj_dir, 'replication/NIH/results');
mshbm_ref_dir = fullfile(proj_dir, 'replication/mshbm');

% Group prior paths
prior_paths = {fullfile(mshbm_ref_dir, 'HCP/Params_Final.mat'), ...
               fullfile(mshbm_ref_dir, 'Du/Params_Final.mat'), ...
               fullfile(mshbm_ref_dir, 'NIH/Params_Final.mat')};

addpath(proj_dir);
addpath(fullfile(proj_dir, 'replication/NIH'));

% MS-HBM scheme parameters
mshbm_flags  = {'HCP', 'Du', 'NIH'};
langA_labels = [5, 9, 8];
langB_labels = [4, 3, 11];
run_pattern  = 'bld([\d]+)';
w = '80';
c = '10';

% ====================================
%% Load Subjects from LI Template
% ====================================

LI_template_path = fullfile(analysis_dir, 'NIH_lang_43_subid.csv');
LI = readtable(LI_template_path);

% Keep only subjects with known language laterality (L, B, or R)
LI = LI(ismember(LI.task, {'L', 'B', 'R'}), :);
fprintf('\nSubjects with valid language task (L/B/R): %d\n', height(LI));

% ====================================
%% Generate Full Individual Parcellations
% ====================================

fprintf('\n====== Generating full individual parcellations ======\n');

for sub_n = 1:height(LI)
    subject    = LI.subid{sub_n};
    surf_files = dir(fullfile(analysis_dir, subject, 'surf', 'lh*sm6.nii.gz'));

    if isempty(surf_files)
        warning('No surf files for %s; skipping.', subject);
        continue
    end

    surf_names  = cellfun(@(f) f(4:end), {surf_files.name}, 'UniformOutput', false);
    lh_list     = cellfun(@(f) fullfile(analysis_dir, subject, 'surf', ['lh.' f]), ...
        surf_names, 'UniformOutput', false);
    rh_list     = cellfun(@(f) fullfile(analysis_dir, subject, 'surf', ['rh.' f]), ...
        surf_names, 'UniformOutput', false);
    censor_list = cellfun(@(f) get_censor_path(analysis_dir, subject, f, run_pattern), ...
        surf_names, 'UniformOutput', false);

    for mshbm_n = 1:numel(mshbm_flags)
        mshbm_flag = mshbm_flags{mshbm_n};
        out_dir    = fullfile(analysis_dir, subject, 'mshbm', mshbm_flag);
        parcel_out = fullfile(out_dir, 'ind_parcellation/test_set', ...
            'Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

        if exist(parcel_out, 'file')
            fprintf('Exists: %s %s, skipping\n', subject, mshbm_flag);
            continue
        end

        params              = struct();
        params.project_dir  = out_dir;
        params.group_prior  = prior_paths{mshbm_n};
        params.lh_fMRI_list = lh_list;
        params.rh_fMRI_list = rh_list;
        params.censor_list  = censor_list;
        params.target_mesh  = 'fsaverage6';
        params.w            = w;
        params.c            = c;

        CBIG_MSHBM_Epilepsy_LI(params);
        fprintf('Done: %s %s\n', subject, mshbm_flag);
    end
end

% ====================================
%% Compute Language Laterality Index
% ====================================

fprintf('\n====== Computing language laterality index ======\n');

for mshbm_n = 1:numel(mshbm_flags)
    mshbm_flag = mshbm_flags{mshbm_n};
    label1     = langA_labels(mshbm_n);
    label2     = langB_labels(mshbm_n);
    col_lang     = [mshbm_flag '_lang'];
    LI.(col_lang) = NaN(height(LI), 1);

    for sub_n = 1:height(LI)
        subject     = LI.subid{sub_n};
        parcel_file = fullfile(analysis_dir, subject, 'mshbm', mshbm_flag, ...
            'ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w80_MRF10.mat');

        if ~exist(parcel_file, 'file')
            continue
        end

        ip     = load(parcel_file);
        lh_lbl = double(ip.lh_labels);
        rh_lbl = double(ip.rh_labels);

        lh_lang = sum(lh_lbl == label1 | lh_lbl == label2);
        rh_lang = sum(rh_lbl == label1 | rh_lbl == label2);

        if (lh_lang + rh_lang) > 0
            LI.(col_lang)(sub_n) = (lh_lang - rh_lang) / (lh_lang + rh_lang);
        end
    end
end

% ====================================
%% Save LI Results
% ====================================

lang_dir = fullfile(results_dir, 'lang');
if ~exist(lang_dir, 'dir')
    mkdir(lang_dir);
end
LI_out_path = fullfile(lang_dir, 'NIH_language_prediction_results.xlsx');
lang_cols_save = strcat(mshbm_flags, '_lang');
LI_out = LI(:, [{'subid', 'task', 'n_runs'}, lang_cols_save]);
writetable(LI_out, LI_out_path);
fprintf('\nLI results saved to: %s\n', LI_out_path);

% ====================================
%% Statistics
% ====================================

% Reload from saved file; remove rows where any lang value is missing
LI_data   = readtable(LI_out_path);
lang_cols = strcat(mshbm_flags, '_lang');
LI_valid  = LI_data(~any(ismissing(LI_data(:, lang_cols)), 2), :);
LI_vals   = table2array(LI_valid(:, lang_cols));
LI_task  = cellstr(categorical(LI_valid.task, {'L', 'B', 'R'}));

fprintf('\nSubjects with complete LI data: %d\n', height(LI_valid));

% One-way ANOVA per scheme
p_anova = zeros(1, numel(mshbm_flags));
for k = 1:numel(mshbm_flags)
    p_anova(k) = anova1(LI_vals(:, k), LI_task, 'off');
end
p_fdr = mafdr(p_anova, 'BHFDR', true);

fprintf('\n--- Language LI One-Way ANOVA (uncorrected) ---\n');
for k = 1:numel(mshbm_flags)
    sig = '';
    if p_anova(k) < 0.05; sig = ' *'; end
    fprintf('%s_lang: p=%.4f%s\n', mshbm_flags{k}, p_anova(k), sig);
end

fprintf('\n--- Language LI One-Way ANOVA (BH-FDR corrected, q<0.05) ---\n');
for k = 1:numel(mshbm_flags)
    sig = '';
    if p_fdr(k) < 0.05; sig = ' *'; end
    fprintf('%s_lang: p_fdr=%.4f%s\n', mshbm_flags{k}, p_fdr(k), sig);
end

% ====================================
%% Plot LI Boxplots
% ====================================

task_cats = categorical(LI_task, {'L', 'B', 'R'});
n_schemes = numel(mshbm_flags);

figure;
for k = 1:n_schemes
    ax = subplot(1, n_schemes, k);
    boxplot(LI_vals(:, k), task_cats, 'Colors', 'k', 'MedianStyle', 'line');
    hold on;

    % Individual data points per group
    task_groups = {'L', 'B', 'R'};
    for t = 1:3
        mask = strcmp(LI_task, task_groups{t});
        scatter(repmat(t, sum(mask), 1), LI_vals(mask, k), ...
            25, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
    end
    hold off;

    ylabel('Language laterality index', 'FontName', 'Arial', 'FontSize', 14);
    title(sprintf('%s  p_{FDR}=%.3g', mshbm_flags{k}, p_fdr(k)), ...
        'FontName', 'Arial', 'FontSize', 14);
    set(ax, 'FontName', 'Arial', 'FontSize', 14, 'Box', 'off', ...
        'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');
end
sgtitle('Language Laterality Index by Task Group', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

% ====================================
%% Plot ROC Curves
% ====================================

% Use NIH prior (best-performing scheme) for ROC
nih_idx = strcmp(mshbm_flags, 'NIH');
scores  = LI_vals(:, nih_idx);
classes = {'L', 'B', 'R'};

% One-vs-rest ROC curves
figure; hold on;
roc_colors = lines(3);
for i = 1:3
    bin_labels = strcmp(LI_task, classes{i});
    if strcmp(classes{i}, 'L')
        use_scores = scores;
    else
        use_scores = -scores;
    end
    [X, Y, ~, AUC] = perfcurve(bin_labels, use_scores, true);
    plot(X, Y, 'LineWidth', 2, 'Color', roc_colors(i, :), ...
        'DisplayName', sprintf('%s vs Rest (AUC=%.2f)', classes{i}, AUC));
end
plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
hold off;
xlabel('False positive rate', 'FontName', 'Arial', 'FontSize', 14);
ylabel('True positive rate', 'FontName', 'Arial', 'FontSize', 14);
legend('Location', 'SouthEast', 'FontName', 'Arial', 'FontSize', 12);
axis square;
set(gca, 'FontName', 'Arial', 'FontSize', 14, 'Box', 'off', ...
    'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');
sgtitle('One-vs-Rest ROC Curves (NIH prior, lang LI)', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

% Binary ROC curve (L vs R; B merged into R)
task_bin   = strrep(LI_task, 'B', 'R');
bin_labels = strcmp(task_bin, 'L');
[X, Y, ~, AUC] = perfcurve(bin_labels, scores, true);

figure; hold on;
plot(X, Y, 'LineWidth', 2, 'Color', 'k', ...
    'DisplayName', sprintf('L vs R (AUC=%.2f)', AUC));
plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
hold off;
xlabel('False positive rate', 'FontName', 'Arial', 'FontSize', 14);
ylabel('True positive rate', 'FontName', 'Arial', 'FontSize', 14);
legend('Location', 'SouthEast', 'FontName', 'Arial', 'FontSize', 12);
axis square;
set(gca, 'FontName', 'Arial', 'FontSize', 14, 'Box', 'off', ...
    'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');
sgtitle('Binary ROC Curve (L vs R, NIH prior, lang LI)', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

% ====================================
%% Cleanup Paths
% ====================================

rmpath(proj_dir);
rmpath(fullfile(proj_dir, 'replication/NIH'));

% ====================================
%% Local Functions
% ====================================

function censor_path = get_censor_path(analysis_dir, subject, surf_file, run_pattern)
% Returns the motion censor file path for a given surface fMRI filename.
bld = regexp(surf_file, run_pattern, 'match');
censor_path = fullfile(analysis_dir, subject, 'qc', ...
    [subject '_' bld{1} '_FDRMS0.2_DVARS50_motion_outliers.txt']);
end
