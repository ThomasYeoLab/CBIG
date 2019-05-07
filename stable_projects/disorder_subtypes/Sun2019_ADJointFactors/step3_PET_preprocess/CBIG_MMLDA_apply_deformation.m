% List of open inputs
% Apply deformations (many subjects): Deformation Field - cfg_files
% Apply deformations (many subjects): Images - cfg_files
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
addpath(spm_dir)
disp(deformation_field)
disp(image_path)

nrun = 1; % enter the number of runs here
jobfile = {[script_dir '/CBIG_MMLDA_apply_deformation_job.m']};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = {deformation_field}; % Apply deformations (many subjects): Deformation Field - cfg_files
    inputs{2, crun} = {image_path}; % Apply deformations (many subjects): Images - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});

rmpath(spm_dir)