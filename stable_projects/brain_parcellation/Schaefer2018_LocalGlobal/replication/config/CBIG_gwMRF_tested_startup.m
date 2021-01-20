% At startup, matlab will find a file startup.m assumed to be located in directory specified by 
% the environmental variable MATLABPATH (see http://www.mathworks.com/help/matlab/ref/startup.html).
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if (~isdeployed)

    if (~isempty(getenv('CBIG_CODE_DIR')))

        % add all paths in CBIG_CODE_DIR subdirectories
        addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab')));

        % add freesurfer/fsfast paths
        addpath(fullfile(getenv('FREESURFER_HOME'), 'fsfast', 'toolbox'));
        addpath(fullfile(getenv('FREESURFER_HOME'), 'matlab'));

        % add external_packages/matlab
        startup_path = pwd;
        cd(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab'));
        external_matlab_add_all_path;
        cd(startup_path);

        % add SD
        CBIG_SD_DIR = getenv('CBIG_SD_DIR');
        if (~isempty(CBIG_SD_DIR) && exist(CBIG_SD_DIR, 'dir'))
            addpath(genpath(CBIG_SD_DIR));
        else
            disp('CBIG_SD_DIR is not set or points to a non-existing directory');
        end
        clear CBIG_SD_DIR;

    else
        disp('CBIG_CODE_DIR not defined!');
    end
    clear all;
end

