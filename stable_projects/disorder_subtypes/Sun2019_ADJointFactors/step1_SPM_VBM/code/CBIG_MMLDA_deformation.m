% List of open inputs
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
addpath(spm_dir)

nrun = 1; % enter the number of runs here
jobfile = {[script_dir '/CBIG_MMLDA_deformation_job.m']};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);
for crun = 1:nrun
	inputs{1, crun} = {[output_dir '/mri/y_Template_6_2mni.nii']};
	inputs{2, crun} = {
	                   [output_dir '/mri/Template_0.nii']
	                   [output_dir '/mri/Template_1.nii']
	                   [output_dir '/mri/Template_2.nii']
	                   [output_dir '/mri/Template_3.nii']
	                   [output_dir '/mri/Template_4.nii']
	                   [output_dir '/mri/Template_5.nii']
	                   [output_dir '/mri/Template_6.nii']
	                   };
	inputs{3, crun} = {[output_dir '/mri']};
	inputs{4, crun} = {[output_dir '/mri/Template_6.nii']};

end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

rmpath(spm_dir)