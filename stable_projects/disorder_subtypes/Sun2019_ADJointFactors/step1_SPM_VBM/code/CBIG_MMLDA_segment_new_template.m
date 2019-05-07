% List of open inputs
% CAT12: Segmentation: Volumes - cfg_files
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
addpath(spm_dir)

nrun = 1; % enter the number of runs here
jobfile = {[script_dir '/CBIG_MMLDA_segment_new_template_job.m']};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);
for crun = 1:nrun
    inputs{1, crun} = {image_path}; % CAT12: Segmentation: Volumes - cfg_files
	inputs{2, crun} = {[spm_dir '/tpm/TPM.nii']};
	inputs{3, crun} = {new_template};
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

rmpath(spm_dir)