function tfcestat = palm_tfce(X,y,opts,plm)
% Compute the TFCE statistic, for volume or surface
% data (vertexwise or facewise).
%
% Usage:
% tfcestat = palm_tfce(X,y,opts,plm)
%
% Inputs:
% - X    : Statistical map.
% - y    : Modality index (of those stored in the plm struct).
% - opts : Struct with PALM options.
% - plm  : Struct with PALM data.
%
% Outputs:
% - tfcestat  : TFCE map.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Sep/2013
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

% Choose an appropriate mask struct.
if opts.npcmod || opts.MV || opts.forcemaskinter,
    S = plm.maskinter;
else
    if plm.nmasks == 1,
        S = plm.masks{1};
    else
        S = plm.masks{y};
    end
end

% Inject the data.
mask    = S.data;
D       = double(S.data);
D(mask) = X;

% Split the code according to whether dh is "automatic" or a fixed
% value given supplied by the user.
if opts.tfce.deltah == 0,
    
    % "delta h"
    dh = max(X(:))/100;
    
    if plm.Yisvol(y),
        
        % Volume (voxelwise data)
        tfcestat = zeros(size(D));
        for h = dh:dh:max(D(:));
            CC    = bwconncomp(D>=h,opts.tfce.conn);
            integ = cellfun(@numel,CC.PixelIdxList).^opts.tfce.E * h^opts.tfce.H;
            for c = 1:CC.NumObjects,
                tfcestat(CC.PixelIdxList{c}) = ...
                    tfcestat(CC.PixelIdxList{c}) + integ(c);
            end
        end
        
    elseif plm.Yisvtx(y) || plm.Yisfac(y),
        
        % Vertexwise or facewise surface data
        tfcestat = zeros(size(D));
        for h = dh:dh:max(D(:));
            dpxl  = palm_dpxlabel(D>=h,plm.Yadjacency{y});
            U     = unique(dpxl(dpxl>0))';
            for u = 1:numel(U),
                idx = dpxl == U(u);
                tfcestat(idx) = tfcestat(idx) + ...
                    sum(plm.Yarea{y}(idx)).^opts.tfce.E * h^opts.tfce.H;
            end
        end
    end
    
else
    if plm.Yisvol(y),
        
        % "delta h"
        dh = opts.tfce.deltah;
        
        % Volume (voxelwise data)
        tfcestat  = zeros(size(D));
        h         = dh;
        CC        = bwconncomp(D>=h,opts.tfce.conn);
        while CC.NumObjects,
            integ = cellfun(@numel,CC.PixelIdxList).^opts.tfce.E * h^opts.tfce.H;
            for c = 1:CC.NumObjects,
                tfcestat(CC.PixelIdxList{c}) = ...
                    tfcestat(CC.PixelIdxList{c}) + integ(c);
            end
            h     = h + opts.tfce.deltah;
            CC    = bwconncomp(D>=h,opts.tfce.conn);
        end
        
    elseif plm.Yisvtx(y) || plm.Yisfac(y),
        
        % Vertexwise or facewise surface data
        tfcestat  = zeros(size(D));
        h         = opts.tfce.deltah;
        dpxl      = palm_dpxlabel(D>=h,plm.Yadjacency{y});
        U         = unique(dpxl(dpxl>0))';
        while numel(U),
            for u = 1:numel(U),
                idx = dpxl == U(u);
                tfcestat(idx) = tfcestat(idx) + ...
                    sum(plm.Yarea{y}(idx)).^opts.tfce.E * h^opts.tfce.H;
            end
            h     = h + opts.tfce.deltah;
            dpxl  = palm_dpxlabel(D>=h,plm.Yadjacency{y});
            U     = unique(dpxl(dpxl>0))';
        end
    end
end

% Return as a vector with the same size as X, and
% apply the correction for the dh.
tfcestat = tfcestat(mask);
tfcestat = tfcestat(:)' * dh;
