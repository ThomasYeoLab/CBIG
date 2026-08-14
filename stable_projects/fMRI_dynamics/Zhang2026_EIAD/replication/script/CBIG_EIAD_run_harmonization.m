function [EI_harmonized, subject_list] = CBIG_EIAD_run_harmonization(scanner_list_dir, cov_file, ...
    S_E_path, S_I_path, r_E_path, combat_dir)
% CBIG_EIAD_RUN_HARMONIZATION  Harmonize regional E/I ratios across scanners with ComBat.
%
% For every participant listed in the scanner lists, this function loads the
% pFIC simulation output, discards the parameter draws whose excitatory firing
% rate falls outside the physiological range, averages the remaining draws to
% obtain the regional excitatory (S_E) and inhibitory (S_I) synaptic gating
% variables, forms the regional E/I ratio (S_E / S_I), and passes the resulting
% 68 x N matrix to ComBat. Scanner identity is the batch variable; age, sex,
% diagnosis, the two age-by-diagnosis interactions and head motion are retained
% as biological covariates so that ComBat removes scanner effects only.
%
% Inputs:
%   scanner_list_dir : folder containing one .txt file per scanner. Each file
%                      lists the subject IDs acquired on that scanner, one per
%                      line. Files are read in the order returned by dir(), and
%                      the k-th file defines batch k.
%   cov_file         : participant table (.xlsx or .csv) whose first five
%                      columns are ID, age, sex (0/1), diagnosis (1 = CN,
%                      2 = MCI, 3 = AD) and mean framewise displacement.
%   S_E_path         : path template for the S_E file of one participant, with a
%                      '%s' placeholder for the subject ID, e.g.
%                      '/my/path/EI/%s/test/simulation/S_E_1.mat'.
%                      The .mat file must contain a 68 x nDraws variable 'S_E'.
%   S_I_path         : same, for the variable 'S_I'.
%   r_E_path         : same, for the excitatory firing rate 'r_E'.
%   combat_dir       : (optional) folder containing combat.m, from
%                      https://github.com/Jfortin1/ComBatHarmonization
%                      (Matlab/scripts). Added to the path when supplied.
%
% Outputs:
%   EI_harmonized : 68 x N matrix of harmonized regional E/I ratios.
%   subject_list  : N x 1 cell array of subject IDs, in the same column order
%                   as EI_harmonized.
%
% Note: the participant-level simulation output (S_E, S_I, r_E) and the
% covariate table contain individual-level data and are not redistributed with
% this repository; this function is provided so that the harmonization
% procedure is fully specified. The harmonized values it produces are the
% 'HarmonizedEI1' ... 'HarmonizedEI68' columns of ADNI_df.csv.
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 6, combat_dir = ''; end
if ~isempty(combat_dir), addpath(combat_dir); end

try
    % firing-rate window used to reject non-physiological parameter draws (Hz)
    r_E_min = 2.7;
    r_E_max = 3.3;

    cov = readtable(cov_file);

    file_list = dir(fullfile(scanner_list_dir, '*.txt'));
    file_list = file_list(~[file_list.isdir]);

    alpha        = [];   % 68 x N regional E/I ratios
    batch        = [];   % 1 x N scanner index
    covariate    = [];   % N x 4 (age, sex, diagnosis, mean FD)
    subject_list = {};

    for i = 1:numel(file_list)
        fprintf('Processing scanner file %d/%d: %s\n', i, numel(file_list), file_list(i).name);

        fid = fopen(fullfile(scanner_list_dir, file_list(i).name), 'r');
        scanner_subjects = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        scanner_subjects = scanner_subjects{1};

        subject_list = [subject_list; scanner_subjects]; %#ok<AGROW>

        for j = 1:numel(scanner_subjects)
            subject = scanner_subjects{j};

            ind = find(strcmp(subject, cov.ID));
            if numel(ind) ~= 1
                error('Subject %s matches %d rows of %s (expected exactly 1).', ...
                    subject, numel(ind), cov_file);
            end

            r_E = load_sim_variable(r_E_path, subject, 'r_E');
            S_E = load_sim_variable(S_E_path, subject, 'S_E');
            S_I = load_sim_variable(S_I_path, subject, 'S_I');

            % drop parameter draws with a non-physiological excitatory firing rate
            anomaly = [find(min(r_E) < r_E_min) find(max(r_E) > r_E_max)];
            S_E(:, anomaly) = [];
            S_I(:, anomaly) = [];

            % average the surviving draws, then form the regional E/I ratio
            SE = mean(S_E, 2, 'omitnan');
            SI = mean(S_I, 2, 'omitnan');

            alpha     = [alpha SE./SI];
            batch     = [batch i];
            covariate = [covariate; cov{ind, 2:5}];
        end
    end

    % Biological covariates retained by ComBat
    age     = covariate(:, 1);
    sex     = dummyvar(covariate(:, 2) + 1);
    sex     = sex(:, 2);                  % 1 = male
    DX      = dummyvar(covariate(:, 3));
    DX      = DX(:, 2:3);                 % MCI and AD indicators (CN is the reference)
    mean_FD = covariate(:, 4);
    age_MCI = age .* DX(:, 1);
    age_AD  = age .* DX(:, 2);

    mod = [age sex DX age_MCI age_AD mean_FD];

    % ComBat harmonization (the trailing 0 selects the non-parametric empirical
    % Bayes adjustment rather than the parametric one)
    EI_harmonized = combat(alpha, batch, mod, 0);

catch ME
    if ~isempty(combat_dir), rmpath(combat_dir); end
    rethrow(ME)
end

if ~isempty(combat_dir), rmpath(combat_dir); end

end

%% ---------------------------------------------------------------------------
function x = load_sim_variable(path_template, subject, var_name)
% Load a single named variable from a per-subject simulation .mat file.
% path_template contains a '%s' placeholder for the subject ID.
f = sprintf(path_template, subject);
if ~exist(f, 'file')
    error('Simulation file not found: %s', f);
end
s = load(f, var_name);
if ~isfield(s, var_name)
    error('Variable ''%s'' not found in %s.', var_name, f);
end
x = s.(var_name);
end
