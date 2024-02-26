function CBIG_preproc_matlab_runtime_compile_utilities(util_file,util_dir,out_dir)

% CBIG_preproc_matlab_runtime_compile_utilities(util_file,util_dir,out_dir)
%
% Compiles MATLAB .m utility function files into standalone executables,
% using the mcc function. The purpose of this is to run these functions
% without installing MATLAB, in MATLAB Runtime. For example, this could be
% useful to run the CBIG preprocessing functions on a computer cluster that
% does not have the paid MATLAB program installed.
%
% Input:
%     - util_file:
%           text file whereby each line contains a single .m utility
%           function, for compilation.
%           For instance, if user intends to compile 2 utility functions:
%           CBIG_preproc_aCompCor.m and CBIG_preproc_censor.m, then
%           util_file should be a text file with two lines. The first line
%           should contain CBIG_preproc_aCompCor.m, and the second line
%           should contain CBIG_preproc_censor.m.
%
%     - util_dir:
%           full path of input directory containing the original .m utility
%           functions that the user intends to compile.
%
%     - out_dir:
%           full path of output directory to contain:
%           1) Compilation log file from running this script for all
%              utility functions. (called "CBIG_matlab_runtime_compile_utilities.log")
%           2) Executable file for each compiled utility function.
%              (for "CBIG_preproc_censor.m", the executable will be named
%              "CBIG_preproc_censor")
%
% Example:
%     util_file = '$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/
%                 matlab_runtime/examples/input/utilities2compile.txt';
%     util_dir = '$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/
%                utilities/'
%     out_dir = '$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/
%                matlab_runtime/utilities/'
%     CBIG_preproc_matlab_runtime_compile_utilities(util_file,util_dir,out_dir)
%
% Written by Trevor Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Read util_file
fid = fopen(util_file);
data = textscan(fid,'%s');
fclose(fid);
util_names = string(data{:}); % create string array of function names

%% Function Compilation
cd(util_dir);
log_file = fullfile(out_dir,'CBIG_matlab_runtime_compile_utilities.log');
mkdir(out_dir);
diary(log_file)
for util_idx = 1:length(util_names)
    util_name = util_names(util_idx);
    mcc('-mv',util_name,'-d',out_dir)
    % -m generates standalone application
    % -v displays verbose output
    % -d specifies output directory
end
diary off

%% Delete Unnecessary Files
cd(out_dir);
delete '*.sh';

end