function palm_ciftiwrite(varargin)
% Provides limited support for writing CIFTI files
% (surface only) using the wb_command as the backend.
%
% C = palm_ciftiwrite(filename,data,extra,wb_command)
%
% Inputs:
% - filename    : Filename of the CIFTI (surface only, dtseries).
% - data        : Actual data that will be saved.
% - extra       : Metadata for the CIFTI.
% - wb_command  : (Optional) Full path to the executable wb_command.
% - toscalar    : (Optional) 
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

% Input arguments
narginchk(1,5);
filename = varargin{1};
data     = varargin{2};
extra    = varargin{3};
if nargin == 4 && ~ isempty(varargin{4}),
    wb_command = varargin{4};
else
    wb_command = 'wb_command';
end
if nargin == 5,
    toscalar = varargin{5};
else
    toscalar = false;
end

% Test the wb_command
try  %#ok
    [status,wb_output] = system(wb_command);
end
if status ~= 0,
    disp(wb_output);
    error('Test failed for your architecture. CIFTI files will not be written.');
end

% Convert to a GIFTI object 
C = gifti([]);
F = fieldnames(extra);
for f = 1:numel(F),
    if ~ strcmpi(F{f},'cifti_file_extension')
        C.private.(F{f}) = extra.(F{f});
    end
end
C.cdata = data;

% Adjust the size in the XML-tree
xt = xmltree(C.private.data{1}.metadata.value);
uid = find(xt,'/CIFTI/Matrix/MatrixIndicesMap');
for u = 1:numel(uid),
    attr = attributes(xt,'get',uid(u));
    for a = 1:numel(attr),
        if strcmp(attr{a}.key,'NumberOfSeriesPoints'),
            xt = attributes(xt,'set',uid(u),a,'NumberOfSeriesPoints',num2str(size(C.cdata,2)));
            break;
        end
    end
end
C.private.data{1}.metadata.value = save(xt);

% Save the GIFTI. This will create a .gii and a .gii.dat.
save(C,strcat(filename,'.gii'),'ExternalFileBinary');

% Include the reset flags to the wb_command (thanks Matt Glasser!)
resetstr = '';
if any(strcmpi(extra.cifti_file_extension,{'dscalar','pscalar'})),
    resetstr = '-reset-scalars';
elseif any(strcmpi(extra.cifti_file_extension,{'dtseries','ptseries'})),
    resetstr = '-reset-timepoints 1 0';
end

% Use the Workbench to convert the GIFTI to CIFTI.
[~] = system(sprintf('%s -cifti-convert -from-gifti-ext %s %s %s',...
    wb_command,strcat(filename,'.gii'),...
    sprintf('%s.%s.nii',filename,extra.cifti_file_extension),resetstr));

% If the CIFTI this is to be converted from timeseries to a scalar
if toscalar && any(strcmpi(extra.cifti_file_extension,{'dtseries','ptseries'})),
    [~] = system(sprintf('%s -cifti-convert-to-scalar %s ROW %s',...
        wb_command,...
        sprintf('%s.%s.nii',filename,extra.cifti_file_extension),...
        sprintf('%s.%s.nii',filename,strrep(extra.cifti_file_extension,'tseries','scalar'))));
    delete(sprintf('%s.%s.nii',filename,extra.cifti_file_extension));
end

% Delete the temporary GIFTI file.
% Note that here the extension is "dat" because it comes from the GIFTI
% I/O. In the writing function, it's "data", because it comes from the
% Workbench.
delete(strcat(filename,'.gii'));
delete(strcat(filename,'.dat'));
