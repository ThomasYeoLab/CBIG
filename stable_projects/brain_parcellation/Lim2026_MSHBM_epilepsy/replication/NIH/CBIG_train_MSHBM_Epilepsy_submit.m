% ====================================
%% Help
% ====================================
%
% CBIG_train_MSHBM_Epilepsy_submit
%
% Submits CBIG_train_MSHBM_Epilepsy.m as a single PBS job to the headnode.
% Run this script from the compiler node. It resolves all environment
% variables on the compiler node (where they are set), writes them into a
% launcher script, and submits that launcher to PBS.
%
% INPUTS:
%   None. Requires CBIG_CODE_DIR to be set in the current environment
%   before starting MATLAB. All other paths are sourced automatically
%   from replication/config/CBIG_epilepsy_config.sh on the compute node.
%
% OUTPUTS:
%   None. Submits a PBS job; logs are written to:
%   $proj_dir/replication/NIH/logs/
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Setup
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir  = fullfile(CBIG_CODE_DIR, 'stable_projects/brain_parcellation/Lim2026_MSHBM_epilepsy');
log_dir   = fullfile(proj_dir, 'replication/NIH/logs');
if ~exist(log_dir, 'dir'); mkdir(log_dir); end

% ====================================
%% Submit PBS Job
% ====================================

train_script  = fullfile(proj_dir, 'replication/NIH/CBIG_train_MSHBM_Epilepsy.m');
config_script = fullfile(proj_dir, 'replication/config/CBIG_epilepsy_config.sh');
log_stem = fullfile(log_dir, 'train_mshbm');
cmd = ['ssh headnode "$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd ' ...
    '\"source ' config_script ' && ' ...
    'matlab -nodesktop -nosplash -r ' ...
    '\\\"run(''' train_script '''); exit;\\\"\" ' ...
    '-walltime 72:00:00 -ncpus 20 -mem 200G ' ...
    '-name ''train_mshbm_epilepsy'' ' ...
    '-joberr ' log_stem '.err ' ...
    '-jobout ' log_stem '.out"'];
fprintf('Submitting command:\n%s\n', cmd);
[st, out] = system(cmd);
if st ~= 0
    fprintf('[ERROR] Job submission failed:\n%s\n', strtrim(out));
else
    fprintf('Job submitted. Logs at:\n');
    fprintf('  %s.out\n', log_stem);
    fprintf('  %s.err\n', log_stem);
end
