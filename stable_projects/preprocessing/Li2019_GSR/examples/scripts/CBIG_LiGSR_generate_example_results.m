function CBIG_LiGSR_generate_example_results

% CBIG_LiGSR_generate_example_results
% 
% This function calls `CBIG_LiGSR_example_Variance_Component.sh` and
% `CBIG_LiGSR_example_KRR.sh` in the same folder to run the example tests.
% 
% Inputs: no input variable needed.
%
% Outputs: this script creates a folder `output` under
% `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples`.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Jingwei Li

exm_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'examples');
outdir = fullfile(exm_dir, 'output');

% create output dir (IMPORTANT)
if(exist(outdir, 'dir'))
    rmdir(outdir, 's')
end
mkdir(outdir);

%% variance component model
system([fullfile(exm_dir, 'scripts', 'CBIG_LiGSR_example_Variance_Component.sh'), ' ', outdir])

%% kernel regression
system([fullfile(exm_dir, 'scripts', 'CBIG_LiGSR_example_KRR.sh'), ' ', outdir])


end

