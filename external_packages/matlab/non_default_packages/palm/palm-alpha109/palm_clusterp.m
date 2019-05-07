function [maxsize,clstat,sizes] = palm_clusterp(X,y,thr,opts,plm,fastpiv,M,psi,res,m,c,o)
% Compute cluster pivotal statistics, for volume or surface
% data (vertexwise or facewise).
% 
% Usage:
% [maxsize,clstat,sizes] = palm_clusterp(X,y,thr,opts,plm)
% 
% Inputs:
% - X    : Statistical map.
% - y    : Modality index (of those stored in the plm struct).
% - thr  : Cluster-forming threshold.
% - opts : Struct with PALM options.
% - plm  : Struct with PALM data.
% 
% Outputs:
% - maxstat : Largest cluster statistic.
% - clstat  : Thresholded map with the cluster sizes (cluster statistic).
% - sizes   : Vector with all cluster sizes.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2017
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

% Inject the data and threshold it.
mask    = S.data;
D       = double(S.data);
D(mask) = X;
Dt      = D >= thr;
I       = double(mask);
I(mask) = 1:size(X,2);

% Compute the sizes and the statistic
sizes = [];
if plm.Yisvol(y),
    
    % Connected components: bwconncomp is slightly faster and
    % less memory intensive than bwlabel
    CC = bwconncomp(Dt);
    
    % A loop here is 4x faster than cellfun
    sizes = zeros(CC.NumObjects,1);
    for u = 1:CC.NumObjects,
        resc = sum(res(:,I(CC.PixelIdxList{u})),2);
        psic = sum(psi(:,I(CC.PixelIdxList{u})),2);
        sizes(u) = fastpiv(M,psic,resc,y,m,c,o,plm);
    end
    
    % Compute the statistic image (this should be for the 1st perm only)
    if nargout > 1,
        clstat = zeros(size(D));
        for u = 1:CC.NumObjects,
            clstat(CC.PixelIdxList{u}) = sizes(u);
        end
        clstat = clstat(mask)';
    end
    
elseif plm.Yisvtx(y) || plm.Yisfac(y),
    
    % Connected components:
    dpxl  = palm_dpxlabel(Dt,plm.Yadjacency{y});
    
    % Compute the cluster stats
    U     = unique(dpxl(dpxl>0))';
    sizes = zeros(size(U));
    for u = 1:numel(U),
        resc = sum(res(:,I(dpxl == U(u))),2);
        psic = sum(psi(:,I(dpxl == U(u))),2);
        sizes(u) = fastpiv(M,psic,resc,y,m,c,o,plm);
    end
    
    % Compute the statistic image (this is normally for the 1st perm only)
    if nargout > 1,
        clstat = zeros(size(D));
        for u = 1:numel(U),
            clstat(dpxl == U(u)) = sizes(u);
        end
        clstat = clstat(mask)';
    end
end

% In fact, only the max matters, because the uncorrected cluster extent
% doesn't make much sense. Make sure the output isn't empty.
if isempty(sizes),
    sizes = 0;
end
maxsize = max(sizes);