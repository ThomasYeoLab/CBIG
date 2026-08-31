% ====================================
%% Help
% ====================================
%
% CBIG_test_MSHBM_esfmri_glm
%
% This script builds the prewhitened GLM design matrices for each subject
% and run in the es_pp dataset and calls FreeSurfer mri_glmfit to estimate
% stimulation-evoked gamma maps. Prewhitening uses per-vertex AR1
% coefficients from ar_matrix_flob.mat. The FLOB 3-basis HRF is convolved
% with the stimulation boxcar to form the design matrix. GLM outputs
% (gamma maps) are consumed by CBIG_test_MSHBM_esfmri_inhomo.m.
% Last tested on 20 May 2026.
%
% PREREQUISITES:
%   None. This script can be run independently of the homo script.
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   Preprocessed es_pp BOLD surface files (*_rest_residc_fs6_sm6_censored.nii.gz) in:
%   $proj_dir/replication/esfmri/input/es_pp/{subject}/surf/
%   Stimulation events TSV files must be present in:
%   $proj_dir/replication/esfmri/input/es_pp/{subject}/ieeg/
%   Motion censor files must be present in:
%   $proj_dir/replication/esfmri/input/es_pp/{subject}/qc/
%   Pre-computed per-vertex AR1 coefficients: replication/esfmri/input/ar_matrix_flob.mat
%   FLOB HRF basis functions: replication/esfmri/input/hrfbasisfns.txt
%   Stimulation contrast file: replication/esfmri/input/stim_contrast_pw
%
% OUTPUTS:
%   Prewhitened regressors: results/glm/{subject}/{subject}_{bld}_WX{1-4}_{hemi}.nii.gz
%   Prewhitened response:   results/glm/{subject}/{subject}_{bld}_WY_{hemi}.nii.gz
%   Global intercept:       results/glm/{subject}/{subject}_{bld}_X0.mat
%   GLM output:             results/glm/{subject}/{subject}_{bld}_glmfit_flob3_{hemi}/
%   Note: GLM outputs contain patient-derived data and are not included in the code release.
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, ...
    'stable_projects/brain_parcellation/Lim2026_MSHBM_Epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');
if isempty(data_dir)
    config_path = fullfile(proj_dir, 'replication/config/CBIG_epilepsy_config.sh');
    [~, extracted] = system(['bash -c "source ' config_path ' && echo $CBIG_EPILEPSY_DATA_DIR"']);
    data_dir = strtrim(extracted);
end
if isempty(data_dir)
    error('CBIG_EPILEPSY_DATA_DIR is not set. Source replication/config/CBIG_epilepsy_config.sh.');
end
input_dir   = fullfile(data_dir, 'esfmri/input');  % read-only input data
results_dir = fullfile(proj_dir, 'replication/esfmri/results');  % writable outputs
espp_dir    = fullfile(input_dir, 'es_pp');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end

% GLM parameters
flob_path   = fullfile(input_dir, 'hrfbasisfns.txt');
c_path      = fullfile(input_dir, 'stim_contrast_pw');
ar_path     = fullfile(input_dir, 'ar_matrix_flob.mat');
run_pattern = 'bld([\d]+)';

% ====================================
%% Load Subjects and Motion QC
% ====================================

subid_file = fullfile(input_dir, 'esfmri_subid.csv');
fid = fopen(subid_file);
fgetl(fid);
subid_raw = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
subid = subid_raw{1};

mc_struct = load(fullfile(espp_dir, 'mc_test.mat'), 'mc_test');
mc_test   = mc_struct.mc_test;

% ====================================
%% GLM: Build Design Matrices and Run mri_glmfit
% ====================================

% Load FLOB 3-basis HRF basis functions (520 rows x 3 columns, sampled at 0.05s)
flob_text   = fileread(flob_path);
flob_dbl    = str2double(strsplit(strtrim(flob_text)));
gamma_basis = reshape(flob_dbl, 3, 520)';

% Load pre-computed per-vertex AR1 coefficients (averaged across all subjects/runs)
arcorr    = load(ar_path);
arcorr_lh = arcorr.ar_mat_lh(:, 2);
arcorr_rh = arcorr.ar_mat_rh(:, 2);

fprintf('\n====== GLM: Building design matrices and running mri_glmfit ======\n');

for sub_n = 1:numel(subid)
    subject = subid{sub_n};
    subj    = strrep(subject, '-', '_');

    if ~isfield(mc_test, subj) || sum(mc_test.(subj).pass) < 1
        continue
    end

    surf_files = dir(fullfile(espp_dir, subject, 'surf', ...
        'lh.*_rest_residc_fs6_sm6_censored.nii.gz'));
    glm_dir = fullfile(results_dir, 'glm', subject);
    if ~exist(glm_dir, 'dir')
        mkdir(glm_dir);
    end

    for run_n = 1:numel(surf_files)
        if mc_test.(subj).pass(run_n) == 0
            continue
        end

        lh_surf = fullfile(surf_files(run_n).folder, surf_files(run_n).name);
        rh_surf = strrep(lh_surf, '/lh.', '/rh.');

        bld_tok = regexp(surf_files(run_n).name, run_pattern, 'match');
        bld_str = bld_tok{1};
        bld_num = str2double(bld_str(4:end));
        run_str = sprintf('%02d', bld_num);

        stim_path   = fullfile(espp_dir, subject, 'ieeg', ...
            [subject '_ses-postop_task-es_run-' run_str '_events.tsv']);
        censor_path = fullfile(espp_dir, subject, 'qc', ...
            [subject '_' bld_str '_FDRMS0.5_DVARS50_motion_outliers.txt']);
        out_dir = fullfile(glm_dir, [subject '_' bld_str '_glmfit_flob3']);

        if exist([out_dir '_lh'], 'dir') && exist([out_dir '_rh'], 'dir')
            fprintf('Skipping (already exists): %s %s\n', subject, bld_str);
            continue
        end

        % Load BOLD data and reshape to (n_vox x T)
        BOLD_lh = MRIread(lh_surf);
        BOLD_rh = MRIread(rh_surf);
        [x_sz, y_sz, z_sz, T] = size(BOLD_lh.vol);
        n_vox  = x_sz * y_sz * z_sz;
        Y_lh   = reshape(BOLD_lh.vol, n_vox, T);
        Y_rh   = reshape(BOLD_rh.vol, n_vox, T);
        Y_lh   = Y_lh(:, 2:end);
        Y_rh   = Y_rh(:, 2:end);

        % Build FLOB design matrix and censor
        nii_hdr = niftiinfo(lh_surf);
        if strcmp(nii_hdr.TimeUnits, 'Second')
            TR = nii_hdr.PixelDimensions(4);
        elseif strcmp(nii_hdr.TimeUnits, 'Millisecond')
            TR = nii_hdr.PixelDimensions(4) / 1000;
        else
            error('Unrecognised TR units for %s run-%s', subject, run_str);
        end
        num_timepoints = nii_hdr.ImageSize(4);
        stim           = table2array(readtable(stim_path, 'FileType', 'text', 'Delimiter', '\t'));
        valid_frames   = logical(load(censor_path));
        stim_tr        = ceil(stim(:, 1) / TR);
        boxcar         = zeros(1, length(valid_frames));
        boxcar(stim_tr) = 1;

        gamma_interp = []; gamma_norm = []; hrf = [];
        hrf_cut = []; hrf_censored = [];
        for col = 1:size(gamma_basis, 2)
            gamma_interp(:, col) = interp1(gamma_basis(:, col), ...
                1:TR/0.05:length(gamma_basis)); %#ok<AGROW>
            gamma_norm(:, col) = ...
                gamma_interp(:, col) / max(max(gamma_interp)); %#ok<AGROW>
            hrf(:, col)          = conv(gamma_norm(:, col), boxcar); %#ok<AGROW>
            hrf_cut(:, col)      = hrf(1:length(valid_frames), col); %#ok<AGROW>
            hrf_censored(:, col) = hrf_cut(valid_frames, col); %#ok<AGROW>
        end
        X = [ones(num_timepoints, 1), hrf_censored];
        X = X(2:end, :);

        % Per-voxel AR1 prewhitening
        X_lh  = zeros(n_vox, T-1, 4);
        X_rh  = zeros(n_vox, T-1, 4);
        WY_lh = zeros(n_vox, T-1);
        WY_rh = zeros(n_vox, T-1);

        for vox_lh = 1:n_vox
            ar1  = arcorr_lh(vox_lh);
            A    = eye(T) + diag(-ar1 * ones(T-1, 1), -1);
            W    = sqrtm(A * A');
            W    = W(2:end, 2:end);
            WY_lh(vox_lh, :)   = W * Y_lh(vox_lh, :)';
            X_lh(vox_lh, :, :) = W * X;
        end

        for vox_rh = 1:n_vox
            ar1  = arcorr_rh(vox_rh);
            A    = eye(T) + diag(-ar1 * ones(T-1, 1), -1);
            W    = sqrtm(A * A');
            W    = W(2:end, 2:end);
            WY_rh(vox_rh, :)   = W * Y_rh(vox_rh, :)';
            X_rh(vox_rh, :, :) = W * X;
        end

        % Save X0 (global intercept for mri_glmfit --X)
        X0_path = fullfile(glm_dir, [subject '_' bld_str '_X0.mat']);
        X0 = ones(T-1, 1); %#ok<NASGU>
        save(X0_path, 'X0', '-v4');

        % Save prewhitened per-voxel regressors (WX1-4) and response (WY) as nifti
        X_base    = fullfile(glm_dir, [subject '_' bld_str '_WX']);
        Y_lh_path = fullfile(glm_dir, [subject '_' bld_str '_WY_lh.nii.gz']);
        Y_rh_path = fullfile(glm_dir, [subject '_' bld_str '_WY_rh.nii.gz']);

        for col = 1:4
            BOLD_lh.vol = reshape(X_lh(:, :, col), x_sz, y_sz, z_sz, T-1);
            BOLD_rh.vol = reshape(X_rh(:, :, col), x_sz, y_sz, z_sz, T-1);
            MRIwrite(BOLD_lh, [X_base num2str(col) '_lh.nii.gz']);
            MRIwrite(BOLD_rh, [X_base num2str(col) '_rh.nii.gz']);
        end
        BOLD_lh.vol = reshape(WY_lh, x_sz, y_sz, z_sz, T-1);
        BOLD_rh.vol = reshape(WY_rh, x_sz, y_sz, z_sz, T-1);
        MRIwrite(BOLD_lh, Y_lh_path);
        MRIwrite(BOLD_rh, Y_rh_path);

        % Run FreeSurfer mri_glmfit with prewhitened per-voxel regressors
        cmd = ['mri_glmfit --y ' Y_lh_path ' --X ' X0_path ...
            ' --pvr ' X_base '1_lh.nii.gz --pvr ' X_base '2_lh.nii.gz' ...
            ' --pvr ' X_base '3_lh.nii.gz --pvr ' X_base '4_lh.nii.gz' ...
            ' --C ' c_path ...
            ' --cortex --nii.gz --surf fsaverage6 lh --glmdir ' out_dir '_lh'];
        system(cmd);
        cmd = ['mri_glmfit --y ' Y_rh_path ' --X ' X0_path ...
            ' --pvr ' X_base '1_rh.nii.gz --pvr ' X_base '2_rh.nii.gz' ...
            ' --pvr ' X_base '3_rh.nii.gz --pvr ' X_base '4_rh.nii.gz' ...
            ' --C ' c_path ...
            ' --cortex --nii.gz --surf fsaverage6 rh --glmdir ' out_dir '_rh'];
        system(cmd);
        fprintf('GLM done: %s %s\n', subject, bld_str);
    end
end
