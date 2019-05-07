% #!/usr/bin/octave -q
function palm_makeniimask(varargin)
% Make a 3D mask from a 4D NIFTI file, removing all NaN,
% Inf, as well as voxels with constant values along the
% fourth dimension. The main feature of this command is
% that it can operate in large files.
% 
% Usage:
% makeniimask 4dfile.nii mask.nii
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Nov/2013
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

% Since this function may be executed outside PALM,
% it's good to have an extra check
if exist('isoctave','file') == 2,
    isoctfun = 'isoctave';
elseif exist('palm_isoctave','file') == 2,
    isoctfun = 'palm_isoctave';
else
    isoctfun = '0';
end

% If this is Octave
if eval(isoctfun),
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % If running as a script, take the input arguments
    cmdname = program_invocation_name();
    if ~ strcmpi(cmdname(end-5:end),'octave'),
        varargin = argv();
        
        % Add the paths for the NIFTI class files (should be a
        % subdirectory). For Matlab, it's assumed it's already
        % in the path (it'd run as part of PALM).
        funpth = fileparts(mfilename('fullpath'));
        addpath(fullfile(funpth,'niftimatlib','matlab'));
        addpath(fullfile(funpth,'niftimatlib','matlab'));
    end
end

% This is probably redundant but fix a bug in an old Matlab version
nargin = numel(varargin);

% Print usage if no inputs are given
if nargin == 0 || strcmp(varargin{1},'-q'),
    fprintf('Make a 3D mask from a 4D NIFTI file, removing all NaN,\n');
    fprintf('Inf, as well as voxels with constant values along the\n');
    fprintf('fourth dimension. The main feature of this command is\n');
    fprintf('that it can operate in large files.\n');
    fprintf('\n');
    fprintf('Usage:\n');
    fprintf('makeniimask 4dfile.nii mask.nii\n');
    fprintf('\n');
    fprintf('_____________________________________\n');
    fprintf('Anderson M. Winkler\n');
    fprintf('FMRIB / Univ. of Oxford\n');
    fprintf('Nov/2013\n');
    fprintf('http://brainder.org\n');
    return;
end

% Usual argument checking
if nargin ~= 2,
    error('Incorrect number of arguments.');
end

% Read the input 4D (NIFTI class)
nii = nifti(varargin{1});

% Init the array for the mask
mskarr = zeros(nii.dat.dim(1:3));

% For strings of voxels, with predefined Y and Z coords
for a = 1:nii.dat.dim(2),
    for b = 1:nii.dat.dim(3),
        
        % Read from the disk
        I = squeeze(nii.dat(:,a,b,:));
        
        % Find NaNs, Infs and constant values
        inan = any(isnan(I),2);
        iinf = any(isinf(I),2);
        icte = sum(diff(I,1,2).^2,2) == 0;
        
        % Put in the mask array
        mskarr(:,a,b) = ~ (inan | iinf | icte);
    end
end

% Prepare to save and save
dat         = file_array;
dat.fname   = varargin{2};
dat.dim     = size(mskarr);
dat.dtype   = 'uint8-le';
dat.offset  = ceil(348/8)*8;
msk         = nifti;
msk.dat     = dat;
msk.mat     = nii.mat;
create(msk);
msk.dat(:,:,:) = mskarr(:,:,:);
