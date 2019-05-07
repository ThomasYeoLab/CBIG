% List of open inputs
% Smooth: Images to Smooth - cfg_files
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
addpath(spm_dir)

nrun = 1; % enter the number of runs here
jobfile = {[script_dir '/CBIG_MMLDA_smooth_job.m']};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = {modulated_image}; % Smooth: Images to Smooth - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

rmpath(spm_dir)