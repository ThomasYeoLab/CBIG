function [Y,maskstruct,Yisvol,Yissrf,Ykindstr,Ytmp] = palm_ready(Yfile,maskstruct,opts,removecte)
% An intermediate function to read input data (-i) and
% voxelwise EVs (-evperdat).
%
% Usage:
% [Y,maskstruct] = palm_ready(Yfile,maskstruct,opts)
%
% Inputs:
% - Yfile      : Filename to be read.
% - maskstruct : Mask struct.
% - opts       : Variable opts.
% - removecte  : Remove constant values?
%
% Outputs:
% - Yfile      : Data, ready to go to Yset or EVset.
% - maskstruct : Mask struct, updated.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Nov/2015
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

Yisvol   = false;
Yissrf   = false;
Ykindstr = '_dat';

% Read a temporary version
Ytmp     = palm_miscread(Yfile,opts.useniiclass,opts.o,opts.precision);

% If this is 4D read with the NIFTI class, it needs a mask now
if strcmp(Ytmp.readwith,'nifticlass') && ndims(Ytmp.data) == 4,
    if isempty(maskstruct),
        % If a mask hasn't been supplied, make one
        tmpmsk = false(Ytmp.extra.dat.dim(1:3));
        for a = 1:Ytmp.extra.dat.dim(2), % y coords
            for b = 1:Ytmp.extra.dat.dim(3), % z coords
                I = squeeze(Ytmp.extra.dat(:,a,b,:));
                inan = any(isnan(I),2);
                iinf = any(isinf(I),2);
                if removecte,
                    icte = sum(diff(I,1,2).^2,2) == 0;
                else
                    icte = false(size(iinf));
                end
                tmpmsk(:,a,b) = ~ (inan | iinf | icte);
            end
        end
        maskstruct = palm_maskstruct(tmpmsk(:)',Ytmp.readwith,Ytmp.extra);
    else
        % If a mask was supplied, check its size
        if any(Ytmp.extra.dat.dim(1:ndims(maskstruct.data)) ~= size(maskstruct.data)),
            error([...
                'The size of the data does not match the size of the mask:\n' ...
                '- Data file %s\n' ...
                '- Mask file %s'],Yfile,maskstruct.filename);
        end
    end
end

% Now deal with the actual data
if ndims(Ytmp.data) == 2, %#ok
    
    % Transpose if that was chosen.
    if opts.transposedata,
        Ytmp.data = Ytmp.data';
    end
    
    % Not all later functions are defined for file_array class,
    % so convert to double
    if strcmp(Ytmp.readwith,'nifticlass'),
        Ytmp.data = double(Ytmp.data);
    end
    
    % This should cover the CSV files and DPX 4D files that
    % were converted to CSV with 'dpx2csv' and then transposed.
    Y = Ytmp.data;
    
elseif ndims(Ytmp.data) == 4,
    
    % Sort out loading for the NIFTI class
    if strcmp(Ytmp.readwith,'nifticlass'),
        tmpmsk = maskstruct.data(:)';
        
        % Read each volume, reshape and apply the mask
        Y = zeros(size(Ytmp.data,4),sum(tmpmsk));
        for n = 1:size(Ytmp.data,4),
            tmp = Ytmp.extra.dat(:,:,:,n);
            tmp = tmp(:)';
            Y(n,:) = tmp(tmpmsk);
        end
    else
        % If not read with the NIFTI class, get all immediately
        Y = palm_conv4to2(Ytmp.data);
    end
end

% Check if the size of data is compatible with size of mask.
% If read with the NIFTI class, this was already taken care of
% and can be skipped.
if ~ strcmp(Ytmp.readwith,'nifticlass'),
    if ~ isempty(maskstruct),
        if size(Y,2) == 1,
            fprintf('Expanding file %s to match the size of the mask.\n',Yfile);
            Y = repmat(Y,[1 numel(maskstruct.data)]);
        elseif size(Y,2) ~= numel(maskstruct.data),
            error([...
                'The size of the data does not match the size of the mask:\n' ...
                '- Data file %s\n' ...
                '- Mask file %s'],Yfile,maskstruct.filename);
        end
    end
end

% Make a mask that removes constant values, Inf and NaN. This will be
% merged with the user-supplied mask, if any, or will be the sole mask
% available to select the datapoints of interest.
if isempty(maskstruct) ...
        && ndims(Ytmp.data) == 4 ...
        && strcmp(Ytmp.readwith,'nifticlass'),
    maskydat = true(1,size(Y,2));
else
    ynan = any(isnan(Y),1);
    yinf = any(isinf(Y),1);
    if removecte,
        ycte = sum(diff(Y,1,1).^2) == 0;
    else
        ycte = false(size(yinf));
    end
    maskydat = ~ (ynan | yinf | ycte);
end

% Now apply the mask created above and the one supplied by the user
% for this modality. If no masks were supplied, create them, except
% for the NIFTI class, which should have been created above
if strcmp(Ytmp.readwith,'nifticlass'),
    maskstruct.data(maskstruct.data) = maskydat(:);
else
    if isempty(maskstruct),
        maskstruct = palm_maskstruct(maskydat,Ytmp.readwith,Ytmp.extra);
    else
        maskydat = maskstruct.data(:) & maskydat(:);
        maskstruct.data = reshape(maskydat,size(maskstruct.data));
    end
end
Y = Y(:,maskydat);

% Prepare a string with a representative name for the kind of data,
% i.e., voxel for volumetric data,
if nargout > 2,
    switch Ytmp.readwith,
        case {'nifticlass','fs_load_nifti','fsl_read_avw',...
                'spm_spm_vol','nii_load_nii'},
            Yisvol   = true;
            Ykindstr = '_vox';
        case {'fs_read_curv','gifti'},
            Yissrf   = true;
            Ykindstr = '_dpv';
        case 'dpxread',
            Yissrf   = true;
            Ykindstr = '_dpx'; % this may be overriden later if a surface file is supplied
        case 'fs_load_mgh',
            if ndims(Ytmp.data) == 4 && ...
                    size(Ytmp.data,2) == 1 && ...
                    size(Ytmp.data,3) == 1,
                Yissrf   = true;
                Ykindstr = '_dpx'; % this may be overriden later if a surface file is supplied
            else
                Yisvol   = true;
                Ykindstr = '_vox';
            end
        otherwise
            Ykindstr = '_dat';
    end
end

% If the Ytmp will also be returned (used in the palm_mediation). Unless
% this is a single vector, there is no need to keep the actual data; it's
% the mask that matters.
if nargout == 6 && numel(maskstruct.data) > 1,
    Ytmp.data = [];
    Yset = [];
end
