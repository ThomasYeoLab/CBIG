function ext = palm_checkprogs
% Test whether some external programs or toolboxes
% are available, namely, FSL, FreeSurfer, SPM and
% Jimmy Shen's NIFTI toolbox.
%
% ext = checkprogs
% 
% 'ext' is a struct containing one field for each
% of these applications, each being containing
% a 0 (false) or 1 (true) depending on whether
% these programs are available or not.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

persistent palm_extern;
if isempty(palm_extern),
        
    % Check for FSL
    palm_extern.fsl = false;
    fsldir = getenv('FSLDIR');
    if ~isempty(fsldir),
        palm_extern.fsl = true;
        addpath(fullfile(fsldir,'etc','matlab'));
        fprintf('Found FSL in %s\n',fsldir);
    end
    
    % Check for FreeSurfer
    palm_extern.fs  = false;
    fshome = getenv('FREESURFER_HOME');
    if ~isempty(fshome),
        palm_extern.fs = true;
        addpath(fullfile(fshome,'matlab'));
        fprintf('Found FreeSurfer in %s\n',fshome);
    end
    
    % Check for the NIFTI toolbox
    palm_extern.nii = false;
    if exist('load_nii') == 2 && exist('save_nii') == 2, %#ok
        palm_extern.nii = true;
    end
    
    % Check for SPM
    palm_extern.spm = false;
    try %#ok
        spm_check_installation('basic');
        palm_extern.spm = true;
        spmpath = fileparts(which('spm'));
        fprintf('Found SPM in %s\n',spmpath);
    end
    
    % Check for the HCP Workbench
    palm_extern.wb_command = false;
    [status,wb_command] = system('which wb_command');
    if status == 0,
        palm_extern.wb_command = true;
        fprintf('Found HCP Workbench executable in %s',wb_command);
    end
    
    % Check the path of PALM and add the path for the NIFTI and GIFTI I/O.
    palm_extern.palmpath = fileparts(mfilename('fullpath'));
    addpath(fullfile(palm_extern.palmpath,'fileio'));
    addpath(fullfile(palm_extern.palmpath,'fileio','extras'));
end
ext = palm_extern;
