function palm_miscwrite(varargin)
% Write various scalar data formats based on the struct X.
% 
% palm_miscwrite(X);
% 
% See 'palm_miscread.m' for a description of the contents of X.
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

X = varargin{1};

switch lower(X.readwith),
    
    case 'textscan',
        
        % Write a generic text file. Note that this doesn't recover
        % the original file (spaces become newlines).
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.txt');
        end
        fid = fopen(X.filename,'w');
        fprintf(fid,'%s\n',X.data{:});
        fclose(fid);
    
    case {'load','csvread'},
        
        % Write a CSV file.
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.csv');
        end
        dlmwrite(X.filename,X.data,'delimiter',',','precision','%0.4f');
        
    case 'vestread',
        
        % Write an FSL "VEST" file
        palm_vestwrite(X.filename,X.data);
        
    case 'msetread',
        
        % Write an MSET (matrix set) file
        palm_msetwrite(X.filename,X.data);
        
    case 'wb_command',
        
        % Write a CIFTI file using the HCP Workbench
        if nargin == 2,
            toscalar = varargin{2};
        else
            toscalar = false;
        end
        siz = size(X.data);
        if siz(1) == 1 && siz(2) > 1,
            X.data = X.data';
        end
        palm_ciftiwrite(X.filename,X.data,X.extra,[],toscalar);
        
    case 'nifticlass',
        
        % Write using the NIFTI class.
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.nii');
        end
        dat = file_array(  ...
            X.filename,    ...
            size(X.data),  ...
            'FLOAT32-LE',  ...
            ceil(348/8)*8);
        nii      = nifti;
        nii.dat  = dat;
        nii.mat  = X.extra.mat;
        if isfield(X.extra,'mat0')
            nii.mat0 = X.extra.mat0;
        else
            nii.mat0 = X.extra.mat;
        end
        create(nii);
        nii.dat(:,:,:) = X.data(:,:,:);
        
    case 'spm_spm_vol',
        
        % Write NIFTI with SPM
        [~,~,fext] = fileparts(X.filename);
        if isempty(fext),
            X.filename = horzcat(X.filename,'.nii');
        end
        X.extra.fname = X.filename;
        X.extra.dt(1) = spm_type('float32'); % for now, save everything as double
        spm_write_vol(X.extra,X.data);
        
    case 'fs_load_nifti',
        
        % Write NIFTI with FreeSurfer
        X.filename = horzcat(X.filename,'.nii.gz');
        X.extra.hdr.vol = X.data;
        X.extra.hdr.datatype = 64; % for now, save everything as double
        X.extra.hdr.bitpix = 64; % for now, save everything as double
        save_nifti(X.extra.hdr,X.filename);
        
    case 'fsl_read_avw',
        
        % Write NIFTI with FSL
        if ~ isfield(X.extra,'vtype'),
            X.extra.vtype = 'd'; % for now, save everything as double
        end
        % The evalc is needed until that empty 'disp' goes away
        try
            [~] = evalc('save_avw(X.data,X.filename,X.extra.vtype,X.extra.voxsize)');
        catch
            save_avw(X.data,X.filename,X.extra.vtype,X.extra.voxsize);
        end
        
    case 'nii_load_nii',
        
        % Write NIFTI with the NIFTI toolbox
        X.filename = horzcat(X.filename,'.nii');
        X.extra.img = X.data;
        save_nii(X.extra,X.filename);
        
    case 'dpxread',
        
        % Write a DPX (DPV or DPF) file
        palm_dpxwrite(X.filename,X.data,X.extra.crd,X.extra.idx);
        
    case 'srfread',
        
        % Write a SRF file
        srfwrite(X.data.vtx,X.data.fac,X.filename);
        
    case 'fs_read_curv',
        
        % Write a FreeSurfer curvature file
        write_curv(X.filename,X.data,X.extra.fnum);
        
    case 'fs_read_surf',
        
        % Write a FreeSurfer surface file
        write_surf(X.filename,X.data.vtx,X.data.fac);
        
    case 'fs_load_mgh',
        
        % Write a FreeSurfer MGH file
        X.filename = horzcat(X.filename,'.mgz');
        save_mgh(X.data,X.filename,X.extra.M,X.extra.mr_parms);
        
    case 'gifti',
        
        % Write a GIFTI file
        X.filename = horzcat(X.filename,'.gii');
        gii = gifti([]);
        F = fieldnames(X.extra);
        for f = 1:numel(F),
            gii.private.(F{f}) = X.extra.(F{f});
        end
        gii.cdata = X.data';
        save(gii,X.filename,gii.private.data{1}.attributes.Encoding);
end
