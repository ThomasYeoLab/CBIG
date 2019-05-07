function X = palm_miscread(filename,useniiclass,tmppath,precision)
% Read various scalar data formats based on the file extension.
%
% X = palm_miscread(filename,useniiclass,tmppath,precision);
%
% filename    : File to be read.
% useniiclass : True/False. Use the NIFTI class when reading
%               NIFTI files. It requires less memory.
% tmppath     : For some file types, indicate whether to store
%               temporary data.
% precision   : Ensure the output data is 'single' or 'double'
%               precision.
%
% X is a struct that contains the fields:
% X.filename  : Contains the name of the file.
% X.readwith  : This tells which program or function was used
%               to read the data. This is useful when saving the
%               data back, to use a compatible function.
% X.data      : Array with the actual data. The size can vary
%               according to what was read.
% X.extra     : Contain extra information, depending on the kind
%               of data that was read and the function or
%               program used for reading.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2013 (first version)
% Dec/2015 (this version)
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

% Check if the file actually exists before doing anything else
if ~ exist(filename,'file')
    error('File not found: %s',filename);
end

% Store the filename, in case there is a need to overwrite this later
X.filename = filename;

% Take the file extension and try to load accordingly
[~,fnam,fext] = fileparts(X.filename);
fext = strdotsplit(strcat(fnam,fext));

% Some formats use external i/o functions that use random numbers. Save
% current state of the random number generator, then restore at the end.
if palm_isoctave,
    state = rand('state'); %#ok
else
    state = rng;
end

switch lower(fext{end}),
    
    case 'txt',
        
        % Read a generic text file
        X.readwith = 'textscan';
        fid = fopen(X.filename);
        X.data = textscan(fid,'%s');
        X.data = X.data{1};
        fclose(fid);
        X.extra = [];
        
    case 'csv',
        
        % Read a CSV file. It has to contain numeric values only.
        % The command 'csvwrite' is a frontend to 'dlmwrite', which calls
        % 'textscan', which on its turn has a limitation of 100k columns.
        % Using 'load' bypass this issue.
        X.readwith = 'load';
        X.data = load(X.filename);
        X.extra = [];
        
    case {'mat','con','fts','grp'},
        
        % Read an FSL "VEST" file.
        X.readwith = 'vestread';
        [X.data,X.extra.PPH] = palm_vestread(X.filename);
        
    case 'mset',
        
        % Set of matrices
        X.readwith = 'mset';
        X.data = palm_msetread(X.filename);
        
    case 'gz',
        
        % Handle (or not) a gzipped NIFTI or CIFTI file.
        if strcmpi(fext{end-1},'nii'),
            
            if any(strcmpi(fext{end-2},{'dscalar','dtseries','dconn','dlabel','ptseries','merge'})),
                
                % Until CIFTI migrates to HDF5, users will have to uncompress manually.
                error('CIFTI files must be uncompressed before they can be read. Use gunzip and try again.');
                
            else
                % Read as NIFTI proper (not CIFTI)
                extern = palm_checkprogs;
                if useniiclass,
                    error([
                        'Reading of gzipped NIFTI files (.nii.gz) is currently disabled\n' ...
                        'If you are sure that your gzipped files, once uncompressed, are not\n' ...
                        'too large to exceed memory limits, you can include the option ''-noniiclass''\n' ...
                        'in the command line. Otherwise, uncompress manually and try again using\n' ...
                        'as input the .nii files instead.\n' ...
                        'File: %s'],X.filename);
                else
                    if extern.fs,       % Read with FreeSurfer
                        X.readwith = 'fs_load_nifti';
                        X.extra.hdr = load_nifti(X.filename);
                        X.data = X.extra.hdr.vol;
                        X.extra.hdr.vol = [];
                    elseif extern.fsl,  % Read with FSL
                        X.readwith = 'fsl_read_avw';
                        [X.data,X.extra.dims,X.extra.voxsize, ...
                            X.extra.bpp,X.extra.endian] = read_avw(X.filename);
                    else
                        error([
                            'Neither FreeSurfer or FSL were found.\n' ...
                            'To use this data you must do one of:\n' ...
                            '- Make sure FreeSurfer is correctly installed and configured,\n' ...
                            '  and that your ''FREESURFER_HOME'' environmental variable is\n' ...
                            '  properly set;\n' ...
                            '- Make sure FSL is correctly installed and configured,\n' ...
                            '  and that your ''FSLDIR'' environmental variable is\n' ...
                            '  properly set;\n' ...
                            '- If you do not have FSL or FS, uncompress the .gz file\n' ...
                            '  and try again.\n' ...
                            'File: %s'],X.filename);
                    end
                end
            end
        else
            error('Unrecognised format with extension %s%s',fext0,fext);
        end
        
    case {'nii','hdr','img'},
        
        % Check for external loaders
        extern = palm_checkprogs;

        % Handle NIFTI and CIFTI files.
        if strcmpi(fext{end},'nii') && ...
                any(strcmpi(fext{end-1},{'dscalar','dtseries','dconn','dlabel','ptseries','merge','pscalar'})),
            
            % Read a CIFTI file.
            if extern.wb_command,
                X.readwith = 'wb_command';
                [X.data,X.extra] = palm_ciftiread(X.filename,tmppath);
                X.data = X.data';
                X.extra.cifti_file_extension = fext{end-1};
            else
                error('Currently cannot read/write CIFTI files without the HCP Workbench.')
            end
            
        else
            % Read a NIFTI file. Note that this will should not
            % be used for ANALYZE.
            if useniiclass,
                X.readwith = 'nifticlass';
                X.extra = nifti(X.filename);
                X.data = X.extra.dat;
            else
                if extern.fs,       % Read with FreeSurfer
                    X.readwith = 'fs_load_nifti';
                    X.extra.hdr = load_nifti(X.filename);
                    X.data = X.extra.hdr.vol;
                    X.extra.hdr.vol = [];
                elseif extern.spm,  % Read with SPM
                    X.readwith = 'spm_spm_vol';
                    X.extra = spm_vol(X.filename);
                    X.data = spm_read_vols(X.extra);
                elseif extern.fsl,  % Read with FSL
                    X.readwith = 'fsl_read_avw';
                    [X.data,X.extra.dims,X.extra.voxsize, ...
                        X.extra.bpp,X.extra.endian] = read_avw(X.filename);
                elseif extern.nii,  % Read with the NIFTI toolbox
                    X.readwith = 'nii_load_nii';
                    X.extra = load_nii(X.filename);
                    X.data = X.extra.img;
                    X.extra.img = [];
                else
                    error([
                        'No FSL, FreeSurfer or SPM were found.\n' ...
                        'To use this data you must have one of these\n' ...
                        'installed and correctly configured.\n' ...
                        'File: %s\n'],X.filename);
                end
            end
        end
        
    case {'dpv','dpf','dpx'},
        
        % Read a DPV/DPF file, in ASCII
        X.readwith = 'dpxread';
        [X.data,X.extra.crd,X.extra.idx] = palm_dpxread(X.filename);
        
    case 'srf',
        
        % Read a SRF file, in ASCII
        X.readwith = 'srfread';
        [X.data.vtx,X.data.fac] = palm_srfread(X.filename);
        
    case {'area','avg_curv','crv','curv', ...
            'h','k','jacobian_white','mid', ...
            'sulc','thickness','volume'},
        
        % Read a FreeSurfer curvature file
        extern = palm_checkprogs;
        if extern.fs,
            X.readwith = 'fs_read_curv';
            [X.data,X.extra.fnum] = read_curv(X.filename);
        else
            error([
                'FreeSurfer was not found. To use this data, make sure\n' ...
                'that FreeSurfer is correctly installed and configured, and\n' ...
                'that your ''FREESURFER_HOME'' environmental variable is\n' ...
                'properly set.\n' ...
                'File: %s\n'],X.filename);
        end
        
    case {'inflated','nofix','orig','pial', ...
            'smoothwm','sphere','reg','white','white_reg'},
        
        % Read a FreeSurfer surface file
        extern = palm_checkprogs;
        if extern.fs,
            X.readwith = 'fs_read_surf';
            [X.data.vtx,X.data.fac] = read_surf(X.filename);
            X.data.fac = X.data.fac + 1;
        else
            error([
                'FreeSurfer was not found. To use this data, make sure\n' ...
                'that FreeSurfer is correctly installed and configured, and\n' ...
                'that your ''FREESURFER_HOME'' environmental variable is\n' ...
                'properly set.\n' ...
                'File: %s\n'],X.filename);
        end
        
    case {'mgh','mgz'},
        
        % Read a FreeSurfer MGH/MGZ file
        extern = palm_checkprogs;
        if extern.fs,
            X.readwith = 'fs_load_mgh';
            [X.data,X.extra.M,X.extra.mr_parms,X.extra.volsz] = load_mgh(X.filename);
        else
            error([
                'FreeSurfer was not found. To use this data, make sure\n' ...
                'that FreeSurfer is correctly installed and configured, and\n' ...
                'that your ''FREESURFER_HOME'' environmental variable is\n' ...
                'properly set.\n' ...
                'File: %s\n'],X.filename);
        end
        
    case 'gii',
        
        % Read a GIFTI file (no mapped file arrays)
        palm_checkprogs; % ensure GIFTI toolbox in the path
        X.readwith = 'gifti';
        gii = gifti(X.filename);
        if isfield(gii,'cdata'),
            X.data = gii.cdata';
        end
        if isfield(gii,'vertices') && isfield(gii,'faces'),
            X.data.vtx  = gii.vertices;
            X.data.fac  = gii.faces;
            if isfield(gii,'mat'),
                vtx = [X.data.vtx ones(size(X.data.vtx,1),1)];
                vtx = vtx * gii.mat;
                X.data.vtx = vtx(:,1:3);
            end
        end
        for d = numel(gii.private.data):-1:1,
            gii.private.data{d}.data = [];
        end
        X.extra = gii.private;
        
    otherwise
        error('File extension %s not known. Data cannot be loaded\n',fext{end});
end

% Restore the state of the random number generator.
if palm_isoctave,
    rand('state',state); %#ok
else
    rng(state);
end

% Enforce a certain precision defined by the user:
if nargin > 3,
    if strcmpi(precision,'double'),
        X.data = double(X.data);
    elseif strcmpi(precision,'single'),
        X.data = single(X.data);
    end
end

% ==============================================================
function spl = strdotsplit(str)
% Split a string at the dots (.).
idx  = find(str == '.');
idxb = [1 idx+1];
idxe = [idx-1 numel(str)];
spl  = cell(numel(idxb),1);
for s = 1:numel(idxb),
    spl{s} = str(idxb(s):idxe(s));
end
