% List of open inputs
% Population to ICBM Registration: Dartel Template - cfg_files
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
addpath(spm_dir)

nrun = 1; % enter the number of runs here
jobfile = {[script_dir '/CBIG_MMLDA_population_to_ICBM_job.m']};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = {[output_dir '/mri/Template_6.nii']}; % Population to ICBM Registration: Dartel Template - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

rmpath(spm_dir)