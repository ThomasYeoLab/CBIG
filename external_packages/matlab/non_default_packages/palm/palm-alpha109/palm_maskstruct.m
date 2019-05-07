function S = palm_maskstruct(mask,readwith,extra)
% Create a struct for a mask, as if it had been read from a file.
% This is useful to save later the data.
%
% Usage:
% M = maskstruct(mask,readwith,extra)
%
% Inputs:
% mask     : A (1 by m) real array.
% readwith : A string telling which function was used to read
%            original data. See 'miscread.m' for help.
% extra    : A struct that varies according to which function
%            was used to read the data.
%
% Usage:
% S        : A struct derived from the 'extra' argument along
%            The mask itself will then be in M.data.
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

% This is common to all cases below.
S.filename = '';
S.readwith = readwith;

% If more information is given, use to create the structs
% that will be used later to mask and save files. Note that
% the field 'filename' should remain empty.
switch lower(readwith),
    
    case {'load','csvread','vestread'},
        
        % If the original data is a CSV or VEST file.
        S.data = mask;
        S.extra = extra;
        
    case 'wb_command',
        S.data  = mask;
        S.extra = extra;
        
    case 'nifticlass',
        
        % If the original data is NIFTI and was read witht the NIFTI class
        S.data          = palm_conv2to4(mask,extra.dat.dim(1:3));
        S.extra.mat     = extra.mat;
        S.extra.dat.dim = extra.dat.dim;
        
    case 'spm_spm_vol',
        
        % If the original data is NIFTI and was read with SPM.
        S.data           = palm_conv2to4(mask,extra(1).dim(1:3));
        S.extra          = extra(1);
        S.extra.dt(1)    = spm_type('float64');
        S.extra.pinfo(1) = 1;
        
    case 'fs_load_nifti',
        
        % If the original data is NIFTI and was read with FreeSurfer.
        S.data                 = palm_conv2to4(mask,extra.hdr.dim(2:4));
        S.extra                = extra;
        S.extra.hdr.scl_slope  = 1;
        S.extra.hdr.dim([1 5]) = [3 1];
        S.extra.hdr.pixdim(5)  = 0;
        S.extra.hdr.datatype   = 64;
        S.extra.hdr.bitpix     = 64;
        
    case 'fsl_read_avw',
        
        % If the original data is NIFTI and was read with FSL.
        S.data  = palm_conv2to4(mask,extra.dims(1:3));
        S.extra = extra;
        if ~ isfield(S.extra,'vtype'),
            S.extra.vtype = 'd';
        end
        
    case 'nii_load_nii',
        
        % If the original data is NIFTI and was read with the NIFTI toolbox.
        S.data                      = palm_conv2to4(mask,extra.hdr.dime.dim(2:4));
        S.extra                     = extra;
        S.extra.hdr.dime.dim([1 5]) = [3 1];
        S.extra.hdr.dime.pixdim(5)  = 0;
        S.extra.hdr.dime.datatype   = 64;
        S.extra.hdr.dime.bitpix     = 64;
        
    case {'fs_read_curv','dpxread'},
        
        % If the original data is an FS curvature.
        S.data  = mask;
        S.extra = extra;
        
    case 'fs_load_mgh',
        
        % If the original data is an FS MGH/MGZ file.
        S.data  = palm_conv2to4(mask,extra.volsz(1:3));
        S.extra = extra;
        
    case 'gifti',
        
        % If the original data is a GIFTI file.
        S.data  = mask;
        S.extra = extra;
        S.extra.data = S.extra.data(1);

    otherwise
        error('Unknown format.')
end
