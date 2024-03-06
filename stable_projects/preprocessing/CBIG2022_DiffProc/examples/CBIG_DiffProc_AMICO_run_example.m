function CBIG_DiffProc_AMICO_run_example(output_dir)
% This function runs the example data for AMICO example.
%
% Inputs:
% -output_dir
%  Directory in which to save results (e.g. fullfile(ref_dir, 'output'))
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% prepare directories
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
    'CBIG2022_DiffProc');
addpath(fullfile(DIFF_CODE_DIR,'examples'));
ref_dir = fullfile(DIFF_CODE_DIR,'examples','AMICO_example_data');

% create the output directory if needed
if exist(output_dir, 'dir')
    rmdir(output_dir, 's')
end
mkdir(output_dir)
                       
% run AMICO - a job will be sent to the scheduler
command = [fullfile(DIFF_CODE_DIR, 'AMICO', 'CBIG_DiffProc_runAMICO.sh'), ' -s ', ...
    fullfile(ref_dir, 'AMICO_subs.txt'), ...
    ' -d ', fullfile(ref_dir, 'diffusion_data'), ' -o ', output_dir, ...
    ' -m ', fullfile(ref_dir, 'b0_mask'), ' -p CBIG_py3'];
[status, cmdout] = system(command);

% remove path
rmpath(fullfile(DIFF_CODE_DIR,'examples'));
end
