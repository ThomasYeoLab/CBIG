function CBIG_DiffProc_TBSS_run_example(output_dir)
% This function runs the example data for TBSS example.
%
% Inputs:
% -output_dir
%  Directory in which to save results (e.g. fullfile(ref_dir, 'output')) 
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% prepare directories 
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
    'CBIG2022_DiffProc');
addpath(fullfile(DIFF_CODE_DIR,'examples'));
ref_dir = fullfile(DIFF_CODE_DIR,'examples','TBSS_example_data');

% create the output directory if needed
if exist(output_dir, 'dir')
    rmdir(output_dir, 's')
end
mkdir(output_dir)

% generate TBSS - this is the same as running the FSL commands in your console
curr_dir = pwd;
copyfile(fullfile(ref_dir,'input'), output_dir);
cd(output_dir)
command = 'tbss_1_preproc *.nii.gz';
[status, cmdout] = system(command);
command = 'tbss_2_reg -T';
[status, cmdout] = system(command);
command = 'tbss_3_postreg -S';
[status, cmdout] = system(command);
command = 'tbss_4_prestats 0.2';
[status, cmdout] = system(command);
command = 'tbss_non_FA MD';
[status, cmdout] = system(command);
cd(curr_dir);

% remove path
rmpath(fullfile(DIFF_CODE_DIR,'examples'));       
end
