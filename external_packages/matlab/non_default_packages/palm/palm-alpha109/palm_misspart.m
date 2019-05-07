function [X,Z,eCm,eCx,Y,imov,ifix,isdiscrete,istwotail] = palm_misspart(M,C,meth,my,mm,mcar,rmethod)
% Partition the model for missing data, generating all the sub-models
% that can later be subjected to NPC.
%
% Usage:
% [X,Z,eCm,eCx,Y,Pidx,isdiscrete] = palm_misspart(M,C,meth,my,mm)
%
% Inputs:
% M          : Design matrix, to be partitioned.
% C          : Contrast that will define the partitioning.
% meth       : Method for the partitioning. It can be any of those available
%              in palm_partition.
% my         : Missing data indicators for the observations (Y).
% mm         : Missing data indicators for the design (M).
% mcar       : Boolean, indicating whether the missing data process is completely
%              at random (true) or not (false).
% rmethod    : Regression & permutation strategy.
%
% Outputs:
% X          : Cell array with sets of EVs of interest.
% Z          : Cell array with sets of nuisance EVs.
% eCm        : Cell array of effective contrasts.
% eCx        : Same as above, but considering only X.
% Y          : Cell array with indices (logical) or data for regression (double).
%              If empty, it's equivalent to a vector index full of ones.
% imov       : Cell array of indices to select the movable observations.
% ifix       : Cell array of indices to select the fixed position observations.
% isdiscrete : Vector indicating if the respective cell array contains both
%              discrete (binary) X and Y, such that the Chi^2 (Yates) test
%              can be be performed.
% istwotail  : Whether the partial test should be run as two-tailed.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2016 Anderson M. Winkler
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

% Initial partitioning of the model:
[x,z,ecm,ecx] = palm_partition(M,C,'guttman');

% Partitioning to obtain the indices:
[mx,mz] = palm_partition(mm,C,'guttman');
if size(C,2) > 1 || any(mx(:) == 2),
    mx = all(mx,2)*2;
    mx2tail = true;
else
    mx2tail = false;
end

% Two-tailedness for missingness in Y:
if any(my == 2),
    my2tail = true;
else
    my2tail = false;
end

% These are (logical) indices for the variables that are available, ie, not missing:
iy = ~ my;
ix = ~ any(mx,2); % if mx is empty, ix is also empty
iz = ~ any(mz,2); % if mz is empty, iz is also empty
ia = true(size(iy));

% These are the actual missing indicators (double) that go in the design:
my = double(my);
mx = double(mx);
mz = double(mz);

% It's simpler and faster to fork the code 16 times than have
% multiple loops and conditions for the up-to-eight equations.
if isempty(mz),
    if       all(iy) &&   all(ix), % Case 1
        idx{1}     = [ia ia ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
    elseif ~ all(iy) &&   all(ix), % Case 2
        if mcar,
            idx{1} = [iy ia ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
        else
            idx{1} = [iy ia ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
            idx{2} = [ia ia ia];   Y{2} = my;   X{2} = x;    Z{2} = [];   eC{2} = ecm;
            istwotail = [false my2tail];
        end
    elseif   all(iy) && ~ all(ix), % Case 3
        if mcar,
            idx{1} = [ia ix ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
        else
            idx{1} = [ia ix ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
            idx{2} = [ia ia ia];   Y{2} = [];   X{2} = mx;   Z{2} = [];   eC{2} = mkcon(X{2},Z{2});
            istwotail = [false mx2tail];
        end
    elseif ~ all(iy) && ~ all(ix), % Case 5
        if mcar,
            idx{1} = [iy ix ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
        else
            idx{1} = [iy ix ia];   Y{1} = [];   X{1} = x;    Z{1} = [];   eC{1} = ecm;
            idx{2} = [ia ix ia];   Y{2} = my;   X{2} = x;    Z{2} = [];   eC{2} = ecm;
            idx{3} = [iy ia ia];   Y{3} = [];   X{3} = mx;   Z{3} = [];   eC{3} = ecm;
            idx{4} = [ia ia ia];   Y{4} = my;   X{4} = mx;   Z{4} = [];   eC{4} = ecm; % discrete
            istwotail = [false my2tail mx2tail (my2tail|mx2tail)];
        end
    end
else
    if       all(iy) &&   all(ix) &&   all(iz), % Case 1
        idx{1}     = [ia ia ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
    elseif ~ all(iy) &&   all(ix) &&   all(iz), % Case 2
        if mcar,
            idx{1} = [iy ia ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [iy ia ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ia ia];   Y{2} = my;   X{2} = x;    Z{2} = z;    eC{2} = ecm;
            istwotail = [false my2tail];
        end
    elseif   all(iy) && ~ all(ix) &&   all(iz), % Case 3
        if mcar,
            idx{1} = [ia ix ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [ia ix ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ia ia];   Y{2} = [];   X{2} = mx;   Z{2} = z;    eC{2} = mkcon(X{2},Z{2});
            istwotail = [false mx2tail];
        end
    elseif   all(iy) &&   all(ix) && ~ all(iz), % Case 4
        if mcar,
            idx{1} = [ia ia iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [ia ia iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ia ia];   Y{2} = [];   X{2} = x;    Z{2} = mz;   eC{2} = mkcon(ecx,Z{2});
            istwotail = [false false];
        end
    elseif ~ all(iy) && ~ all(ix) &&   all(iz), % Case 5
        if mcar,
            idx{1} = [iy ix ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [iy ix ia];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ix ia];   Y{2} = my;   X{2} = x;    Z{2} = z;    eC{2} = ecm;
            idx{3} = [iy ia ia];   Y{3} = [];   X{3} = mx;   Z{3} = z;    eC{3} = mkcon(X{3},Z{3});
            idx{4} = [ia ia ia];   Y{4} = my;   X{4} = mx;   Z{4} = z;    eC{4} = mkcon(X{4},Z{4}); % discrete
            istwotail = [false my2tail mx2tail (my2tail|mx2tail)];
        end
    elseif ~ all(iy) &&   all(ix) && ~ all(iz), % Case 6
        if mcar,
            idx{1} = [iy ia iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [iy ia iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ia iz];   Y{2} = my;   X{2} = x;    Z{2} = z;    eC{2} = ecm;
            idx{3} = [iy ia ia];   Y{3} = [];   X{3} = x;    Z{3} = mz;   eC{3} = mkcon(ecx,Z{3});
            idx{4} = [ia ia ia];   Y{4} = my;   X{4} = x;    Z{4} = mz;   eC{4} = mkcon(ecx,Z{4});
            istwotail = [false my2tail false my2tail];
        end
    elseif   all(iy) && ~ all(ix) && ~ all(iz), % Case 7
        if mcar,
            idx{1} = [ia ix iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [ia ix iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ia iz];   Y{2} = [];   X{2} = mx;   Z{2} = z;    eC{2} = mkcon(X{2},Z{2});
            idx{3} = [ia ix ia];   Y{3} = [];   X{3} = x;    Z{3} = mz;   eC{3} = mkcon(ecx, Z{3});
            idx{4} = [ia ia ia];   Y{4} = [];   X{4} = mx;   Z{4} = mz;   eC{4} = mkcon(X{4},Z{4});
            istwotail = [false mx2tail false mx2tail];
        end
    elseif ~ all(iy) && ~ all(ix) && ~ all(iz), % Case 8
        if mcar,
            idx{1} = [iy ix iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
        else
            idx{1} = [iy ix iz];   Y{1} = [];   X{1} = x;    Z{1} = z;    eC{1} = ecm;
            idx{2} = [ia ix iz];   Y{2} = my;   X{2} = x;    Z{2} = z;    eC{2} = ecm;
            idx{3} = [iy ia iz];   Y{3} = [];   X{3} = mx;   Z{3} = z;    eC{3} = mkcon(X{3},Z{3});
            idx{4} = [iy ix ia];   Y{4} = [];   X{4} = x;    Z{4} = mz;   eC{4} = mkcon(ecx, Z{4});
            idx{5} = [ia ia iz];   Y{5} = my;   X{5} = mx;   Z{5} = z;    eC{5} = mkcon(X{5},Z{5}); % discrete?
            idx{6} = [ia ix ia];   Y{6} = my;   X{6} = x;    Z{6} = mz;   eC{6} = mkcon(ecx, Z{6});
            idx{7} = [iy ia ia];   Y{7} = [];   X{7} = mx;   Z{7} = mz;   eC{7} = mkcon(X{7},Z{7});
            idx{8} = [ia ia ia];   Y{8} = my;   X{8} = mx;   Z{8} = mz;   eC{8} = mkcon(X{8},Z{8}); % discrete
            istwotail = [false my2tail mx2tail false (my2tail|mx2tail) my2tail mx2tail (my2tail|mx2tail)];
        end
    end
end
nO = numel(idx);

% Overwrite those effective contrast stuff
for c = 1:numel(eC),
    eC{c} = ecm;
end

% For the cases in which there's nothing to be combined, not two-tailed
if nO == 1,
    istwotail = false;
end

% Partition again each of these models using the method indicated by the user:
eCm = cell(size(X));
eCx = eCm;
isdiscrete = false(1,numel(Z));
for o = 1:nO,
    [X{o},Z{o},eCm{o},eCx{o}] = palm_partition(horzcat(X{o},Z{o}),eC{o},meth);
    % Check if discrete. If yes, it will done via Yates' Chi^2.
    % Otherwise, add an intercept
    if ~ isempty(Y{o}) && isempty(Z{o}) && size(unique(X{o},'rows'),1) == 2,
        isdiscrete(o) = true;
    else
        Z{o}   = horzcat(Z{o},ones(size(X{o},1),1));
        eCm{o} = vertcat(eCm{o},zeros(1,size(eCm{o},2)));
    end
end

% Remove bits that are all zeroes (this needs to be adapted for voxelwise):
for o = nO:-1:1,
    if  (~isempty(Y{o}) && ~islogical(Y{o}) && all(Y{o} == 0,1)) || all(X{o}(:,1,1) == 0,1),
        idx(o) = []; Y(o) = []; X(o) = []; Z(o) = [];
        eCm(o) = []; eCx(o) = [];
        isdiscrete(o) = [];
    end
end
nO = numel(isdiscrete);

% Prepare the indices used to modify the permutation matrix
% and model. That is:
imov = cell(nO,1);
ifix = imov;
switch lower(rmethod),
    case 'draper-stoneman',
        for o = 1:nO,
            imov{o} = idx{o}(:,2);
            ifix{o} = all(idx{o}(:,[1 3]),2);
        end
    case 'still-white',
        for o = 1:nO,
            imov{o} = [];
            ifix{o} = all(idx{o}(:,1:2),2);
        end
    case 'freedman-lane',
        for o = 1:nO,
            imov{o} = [];
            ifix{o} = all(idx{o}(:,1:2),2);
        end
    case 'manly',
        for o = 1:nO,
            imov{o} = idx{o}(:,1);
            ifix{o} = all(idx{o}(:,2:3),2);
        end
    case 'terbraak',
        for o = 1:nO,
            imov{o} = [];
            ifix{o} = idx{o}(:,1);
        end
    case 'kennedy',
        for o = 1:nO,
            imov{o} = [];
            ifix{o} = all(idx{o}(:,1:2),2);
        end
    case 'dekker',
        for o = 1:nO,
            imov{o} = [];
            ifix{o} = all(idx{o}(:,1:2),2);
        end
end

% Remove all the full true indices, for speed later:
for o = 1:nO,
    if all(imov{o}),
        imov{o} = [];
    end
    if all(ifix{o}),
        ifix{o} = [];
    end
end

% ==============================================================
function eC = mkcon(X,Z)
% Shortcut to create the effective contrast.
if isempty(Z) || size(X,1) == size(Z,1),
    % Typical case, X is X and Z is Z.
    eC = vertcat(eye(size(X,2)),zeros(size(Z,2),size(X,2)));
else
    % Here X is C, and Z is Z.
    eC = vertcat(X,zeros(size(Z,2),size(X,2)));
end

% ==============================================================
function [Z,eC] = pcaz(Z,eC)
% PCA of Z, via SVD (currently unused).
for o = 1:numel(Z),
    sZ1     = size(Z{o},2);
    Z0      = bsxfun(@minus,Z{o},mean(Z{o},1));
    [u,s,~] = svd(Z0,'econ');
    ds      = diag(s);
    tol     = 100 * max(size(Z{o})) * eps(max(ds));
    abv     = ds > tol;
    Z{o}    = u(:,abv)*s(abv,abv);
    sZ2     = size(Z{o},2);
    eC{o}   = eC{o}(1:end-(sZ1-sZ2),:);
end

% ==============================================================
function [Y,X,Z] = meancenter(Y,X,Z)
% Mean center data and design (currently unused).
for o = 1:numel(Y),
    if ~isempty(Y{o}) && ~islogical(Y{o}),
        Y{o} = bsxfun(@minus,Y{o},mean(Y{o},1));
    end
    X{o} = bsxfun(@minus,X{o},mean(X{o},1));
    Z{o} = bsxfun(@minus,Z{o},mean(Z{o},1));
end