function [data,extra] = palm_ciftiread(varargin)
% Provides limited support for reading CIFTI files
% (surface only) using the wb_command as the backend.
%
% [data,extra] = palm_ciftiread(filename,temp_prefix,wb_command)
%
% Inputs:
% - filename    : Filename of the CIFTI (surface only, dtscalar).
% - temp_prefix : (Optional) Prefix, possibly prepended by a full path
%                 where temporary files will be safed. If omitted,
%                 the current directory will be used.
% - wb_command  : (Optional) Full path to the executable wb_command.
%
% Outputs:
% - data        : Array with the actual data.
% - extra       : Extra information needed to save back as CIFTI.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Apr/2015
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

% Handle input arguments
narginchk(1,3);
filename        = varargin{1};
if nargin >= 2,
    temp_prefix = varargin{2};
else
    temp_prefix = '';
end
if nargin == 3,
    wb_command  = varargin{3};
else
    wb_command  = 'wb_command';
end

% Test the wb_command
try  %#ok
    [status,wb_output] = system(wb_command);
end
if status ~= 0,
    disp(wb_output);
    error('Test failed for your architecture. CIFTI will not be read.');
end

% Deal with the location of the temporary files
alphabet    = ['a':'z' 'A':'Z' '0':'9'];
rand_prefix = alphabet(randi(numel(alphabet),[1 6]));

% If no temporary prefix has been supplied, the file will be
% saved to the current directory, with temp_prefix only.
if ~ isempty(temp_prefix),
    rand_prefix = strcat(temp_prefix,'_',rand_prefix);
end
temp_gii = strcat(rand_prefix,'.gii');

% Convert CIFTI to GIFTI using the Workbench. This line causes troubles in
% Matlab for Mac, because of the paths of the linked libraries of
% the Workbench. Hopefully in the future there will be stable CIFTI I/O
% that will work natively in Matlab and Octave.
% Also note that the wb_command must be accessible with the environmental var PATH.
[~] = system(sprintf('%s -cifti-convert -to-gifti-ext %s %s',wb_command,filename,temp_gii));

% Load the temporary GIFTI file. Note that there are no mapped file-arrays.
gii = gifti(temp_gii);
if isfield(gii,'cdata'),
    data = gii.cdata';
end
if isfield(gii,'vertices') && isfield(gii,'faces'),
    data.vtx = gii.vertices;
    data.fac = gii.faces;
    if isfield(gii,'mat'),
        vtx = [data.vtx ones(size(data.vtx,1),1)];
        vtx = vtx * gii.mat;
        data.vtx = vtx(:,1:3);
    end
end
for d = numel(gii.private.data):-1:1,
    gii.private.data{d}.data = [];
end
extra = gii.private;

% Delete the temporary GIFTI file.
% Note that here the extension is "data" because it comes from the
% Workbench. In the writing function, it's just "dat", because it
% comes from the GIFTI I/O.
delete(temp_gii);
delete(strcat(temp_gii,'.data'));
