function CBIG_TRBPC_example_wrapper(outdir)

% CBIG_TRBPC_example_wrapper(outdir)
%
% This function will run examples sequentially
% 
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


example_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes','ChenTam2022_TRBPC','examples');

% run single-kernel KRR
cmd = [fullfile(example_dir, 'scripts', 'CBIG_TRBPC_example_singleKRR.sh'), ' ', outdir];
system(cmd);

% run multi-kernel KRR
cmd = [fullfile(example_dir, 'scripts', 'CBIG_TRBPC_example_multiKRR.sh'), ' ', outdir];
system(cmd);

% run LRR
cmd = [fullfile(example_dir, 'scripts', 'CBIG_TRBPC_example_LRR.sh'), ' ', outdir];
system(cmd);

% run PFM
cmd = [fullfile(example_dir, 'scripts', 'CBIG_TRBPC_example_PFM.sh'), ' ', outdir];
system(cmd);

% wait for jobs to be finished
cmdout = 1;
while(cmdout ~= 0)
    cmd = 'qstat | grep compute_PFM | grep `whoami` | wc -l';
    [~, cmdout] = system(cmd);
    % after job finishes, cmdout should be 0
    cmdout = str2num(cmdout(1: end-1));
    pause(60); % sleep for 1min and check again
end

end