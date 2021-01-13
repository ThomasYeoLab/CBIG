function [ cifti ] = ciftiopen(filename,wbcommand,verbose)

% CIFTI = CIFTIOPEN(FILENAME,WBCOMMAND,VERBOSE)
%
% Open a CIFTI file by converting to GIFTI external binary first and then
% using the GIFTI toolbox
%
% FILENAME is string containing file name to open
% WBCOMMAND is string containing the Workbench command.  
%   (Necessary for the intermediate step of conversion to gifti.
%    Matlab must be able to find this when it executes a 'system' command).
% VERBOSE (optional; default is off): Set to 1 for more verbose output.

% MPH: Modified to use a temporary file location suitable
% for that OS; added an explicit verbose flag; switched to
% using matlabs 'delete' command for file removal and 'system'
% command for executing 'wbcommand'.

% Default is VERBOSE=0 (OFF)
if (nargin < 3) 
  verbose = 0;
end

% Use Matlab 'tempname' to return a temporary file name and location appropriate for this system
tmpfile = tempname;

tstart=tic;

% Do conversion and loading
system([wbcommand ' -cifti-convert -to-gifti-ext ' filename ' ' tmpfile '.gii']);
cifti = gifti([tmpfile '.gii']);

if (verbose)
  fprintf(1,'%s: Elapsed time is %.2f seconds\n',filename,toc(tstart));
end

% Clean-up
delete([tmpfile '.gii'],[tmpfile '.gii.data']);
