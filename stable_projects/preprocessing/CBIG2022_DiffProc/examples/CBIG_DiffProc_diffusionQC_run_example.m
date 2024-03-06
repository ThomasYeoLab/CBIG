function CBIG_DiffProc_diffusionQC_run_example(output_dir)
% This function runs the example data for diffusionQC example.
%
% Inputs:
% -output_dir
%  Directory in which to save results (e.g. fullfile(ref_dir, 'output'))
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% prepare directories
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); 
DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects','preprocessing', ...
    'CBIG2022_DiffProc');
addpath(fullfile(DIFF_CODE_DIR,'examples'));
ref_dir = fullfile(DIFF_CODE_DIR,'examples','diffusionQC_example_data');
            
% create the output directory if needed
if exist(output_dir, 'dir')
    rmdir(output_dir, 's')
end
            
% run QC pipeline - either run the command directly in your console or
% submit a job to the scheduler
command = strcat(fullfile(DIFF_CODE_DIR, 'CBIG_DiffProc_diffusionQC.sh'), ...
" -s sub-A00059845 -d ", fullfile(ref_dir, 'input', 'sub-A00059845', 'ses-DS2', 'dwi'), ...
" -i sub-A00059845_ses-DS2_dwi.nii.gz -a sub-A00059845_ses-DS2_dwi.bval", ...
" -e sub-A00059845_ses-DS2_dwi.bvec --round_bvals_100 ", ...
" -o ", output_dir);
submit_command = strcat(fullfile(CBIG_CODE_DIR, 'setup', 'CBIG_pbsubmit'), " -cmd '", ...
command, "' -walltime '0:30:00' -name 'DiffQC_example' -mem '8G'");
[status, cmdout] = system(submit_command);

% remove path      
rmpath(fullfile(DIFF_CODE_DIR,'examples'));
end
