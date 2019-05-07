function [Bset,nB,mtr] = palm_shuftree(varargin)
% This is a wrapper for the palm_permtree.m and palm_fliptree.m
% that generates a sigle set of permutations. It can also generate
% only permutations with sign-flipping depending on the input
% arguments.
%
% Usage (style 1)
% [Bset,nB] = palm_shuftree(Ptree,nP0,CMC,EE,ISE,idxout)
%
% Inputs:
% - Ptree   : Permutation tree.
% - nP0     : Requested number of permutations.
% - CMC     : Use Conditional Monte Carlo.
% - EE      : Allow permutations?
% - ISE     : Allow sign-flips?
%             If you supply the EE argument, you must
%             also supply ISE argument. If one is omited,
%             the other needs to be omited too.
%             Default is true for EE, and false for ISE.
% - idxout  : (Optional) If true, the output isn't a cell
%             array with permutation matrices, but an array
%             with permutation indices.
%
% Outputs:
% - Bset    : Set of permutations and/or sign flips.
% - nB      : Number of permutations and/or sign-flips.
%
%
% Usage (style 2, to be used by the PALM main function):
% [Bset,nB,metr] = palm_shuftree(opts,plm)
%
% Inputs:
% - opts    : Struct with PALM options
% - plm     : Struct with PALM data
%
% Outputs:
% - Bset    : Set of permutations and/or sign flips.
% - nB      : Number of permutations and/or sign-flips.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
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

% Take arguments
if nargin == 2 || nargin == 4,
    opts     = varargin{1};
    plm      = varargin{2};
    EE       = opts.EE;
    ISE      = opts.ISE;
    nP0      = opts.nP0;
    CMC      = opts.cmcp;
    seq      = plm.seq{varargin{3}}{varargin{4}};
    Ptree    = palm_tree(plm.EB,seq);
    idxout   = false;
elseif nargin == 3 || nargin == 5 || nargin == 6,
    Ptree    = varargin{1};
    nP0      = varargin{2};
    CMC      = varargin{3};
    if nargin == 5 || nargin == 6,
        EE   = varargin{4};
        ISE  = varargin{5};
    else
        EE   = true;
        ISE  = false;
    end
    if nargin == 6,
        idxout = varargin{6};
    else
        idxout = false;
    end
else
    error('Incorrect number of input arguments');
end
if ~EE && ~ISE,
    error('EE and/or ISE must be enabled, otherwise there is nothing to shuffle.')
end

% Maximum number of shufflings (perms, sign-flips or both)
maxP = 1;
maxS = 1;
if EE,
    lmaxP = palm_maxshuf(Ptree,'perms',true);
    maxP = exp(lmaxP);
    if isinf(maxP),
        fprintf('Number of possible permutations is exp(%g).\n',lmaxP);
    else
        fprintf('Number of possible permutations is %g.\n',maxP);
    end
end
if ISE,
    lmaxS = palm_maxshuf(Ptree,'flips',true);
    maxS = exp(lmaxS);
    if isinf(maxS),
        fprintf('Number of possible sign-flips is exp(%g).\n',lmaxS);
    else
        fprintf('Number of possible sign-flips is %g.\n',maxS);
    end
end
maxB = maxP * maxS;

% String for the screen output below
if EE && ~ISE,
    whatshuf  = 'permutations only';
    whatshuf2 = 'perms';
elseif ISE && ~EE,
    whatshuf  = 'sign-flips only';
    whatshuf2 = 'flips';
elseif EE && ISE,
    whatshuf  = 'permutations and sign-flips';
    whatshuf2 = 'both';
end

% Generate the Pset and Sset
Pset = {};
Sset = {};
if nP0 == 0 || nP0 >= maxB,
    % Run exhaustively if the user requests too many permutations.
    % Note that here CMC is irrelevant.
    fprintf('Generating %g shufflings (%s).\n',maxB,whatshuf);
    if EE,
        Pset = palm_permtree(Ptree,round(maxP),[],false,round(maxP));
    end
    if ISE,
        Sset = palm_fliptree(Ptree,round(maxS),[],false,round(maxS));
    end
elseif nP0 < maxB,
    % Or use a subset of possible permutations. The nested conditions
    % are to avoid repetitions, and to compensate fewer flips with more
    % perms or vice versa as needed in the tight situations
    fprintf('Generating %g shufflings (%s).\n',nP0,whatshuf);
    if EE,
        if nP0 >= maxP,
            Pset = palm_permtree(Ptree,round(maxP),CMC,false,round(maxP));
        else
            Pset = palm_permtree(Ptree,nP0,CMC,false,round(maxP));
        end
    end
    if ISE,
        if nP0 >= maxS,
            Sset = palm_fliptree(Ptree,round(maxS),CMC,false,round(maxS));
        else
            Sset = palm_fliptree(Ptree,nP0,CMC,false,round(maxS));
        end
    end
end

% This ensures that there is at least 1 permutation (no permutation)
% and 1 sign-flipping (no sign-flipping).
nP = numel(Pset);
nS = numel(Sset);
if nP > 0 && nS == 0,
    Sset{1} = Pset{1};
    nS = 1;
elseif nP == 0 && nS > 0,
    Pset{1} = Sset{1};
    nP = 1;
end

% Generate the set of shufflings, mixing permutations and
% sign-flippings as needed.
if nS == 1,
    % If only 1 sign-flip is possible, ignore it.
    Bset = Pset;
elseif nP == 1,
    % If only 1 permutation is possible, ignore it.
    Bset = Sset;
elseif nP0 == 0 || nP0 >= maxB,
    % If the user requested too many shufflings, do all
    % those that are possible.
    Bset = cell(maxB,1);
    b = 1;
    for p = 1:numel(Pset),
        for s = 1:numel(Sset),
            Bset{b} = Pset{p} * Sset{s};
            b = b + 1;
        end
    end
else
    % The typical case, with an enormous number of possible
    % shufflings, and the user choses a moderate number
    Bset = cell(nP0,1);
    % 1st shuffling is no shuffling, regardless
    Bset{1} = Pset{1} * Sset{1};
    if CMC,
        % If CMC, no need to take care of repetitions.
        for b = 2:nP0,
            Bset{b} = Pset{randi(nP)} * Sset{randi(nS)};
        end
    else
        % Otherwise, avoid them
        [~,idx] = sort(rand(nP*nS,1));
        idx = idx(1:nP0);
        [pidx,sidx] = ind2sub([nP nS],idx);
        for b = 2:nP0,
            Bset{b} = Pset{pidx(b)} * Sset{sidx(b)};
        end
    end
end
nB = numel(Bset);

% In the draft mode, the permutations can't be in lexicographic
% order, but entirely shuffled.
if nargin == 2 || nargin == 4,
    if opts.accel.negbin,
        Bset2 = cell(size(Bset));
        [~,idx] = sort(rand(nB,1));
        for p = 2:nB,
            Bset2{p} = Bset(idx(p));
        end
        Bset = Bset2;
    end
end

% If the desired outputs are permutation indices instead
% of permutation matrices, output them
if idxout || ... % indices out instead of a cell array
        (nargout == 3 && nargin == 4), % save metrics
    
    % Convert formats
    Bidx = palm_swapfmt(Bset);
        
    % Compute some metrics
    if nargout == 3,
        Ptree1 = palm_tree(plm.EB,ones(size(seq)));
        mtr = zeros(9,1);
        [mtr(1),mtr(2),mtr(4)] = ...
            palm_metrics(Ptree,seq,whatshuf2);
        [~,~,mtr(3)] = ...
            palm_metrics(Ptree1,ones(size(seq)),whatshuf2);
        [mtr(5),mtr(6),mtr(7),mtr(8),mtr(9)] = palm_metrics(Bidx,seq,whatshuf2);
    end
    
    % Output as indices if needed
    if idxout,
        Bset = Bidx;
    end
end
