function varargout = palm_metrics(varargin)
% Compute some permutation metrics:
% - For permutation trees, return the entropies.
% - For sets of permutations, return the average Hamming distance.
% 
% Usage:
% [lW,lW0,C] = palm_metrics(Ptree,X,stype)
% [Hamm,HammX,Eucl,EuclX,Spear] = palm_metrics(Pset,X)
% 
% Inputs:
% - Ptree : Permutation tree.
% - X     : Design matrix (only the EVs of interest for the Freedman-Lane
%           and most methods, or the full matrix for ter Braak).
%           Note that the metrics are only meaningful if X is the same
%           used when Ptree was created originally.
% - stype : Shuffling type. It can be 'perms', 'flips' or 'both'.
% - Pset  : Set of shufflings (permutations or sign-flips).
% 
% Outputs
% - lW    : Log of the max number of permutations with the restrictions
%           imposed by the tree and the original design used to create the
%           tree.
% - lW0   : Log of the max number of permutations without the restrictions
%           imposed by the tree, but with the restrictions imposed by the
%           input design X.
% - C     : Huberman & Hogg complexity (C) of a given tree.
%           For this to give exactly the same result as in the original
%           paper, such that it measures the relationships in the tree
%           itself, rather than the actual values found in X, the input
%           Ptree must have been constructed with X = ones(N,1) (or any
%           other constant). However, C doesn't depend on the X that is
%           input (i.e., it's not an argument needed to compute C, but
%           it's implicitly taken into account through the tree).
% - Hamm  : Average Hamming distance across the given permutation set,
%           i.e., it's the average change that a permutation cause on
%           the original indices.
% - HammX : Same as Hamm, but consider the possibility of repeated
%           elements in X. If X isn't supplied, or if X has no ties,
%           or if X is the same used originally to create the permutation
%           set, Hamm and HammX are the same.
% - Eucl  : Same as Hamm, but using the Euclidean distance.
% - EuclX : Same as HammX, but using the Euclidean distance.
% - Spear : Same as Hamm, but using the Spearman correlation. X is always
%           taken into account.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2014 (updated Oct/2014)
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

% Take args and decide what to do
X = [];
if iscell(varargin{1}),
    dowhat = 'entropy';
    Ptree  = varargin{1};
    N      = numel(palm_permtree(Ptree,1,false,true));
else
    dowhat = 'distances';
    Pset   = varargin{1};
    N      = size(Pset,1);
end
if nargin == 1,
    X     = (1:N)';
    stype = 'perms';
elseif nargin == 2,
    X     = varargin{2};
    stype = 'perms';
elseif nargin == 3,
    X     = varargin{2};
    stype = varargin{3};
end
if isempty(X),
    X = (1:N)';
end

switch dowhat,
    case 'entropy',
        
        % Normalised entropy or anisotropy. This is computed
        % if the first input is a cell array (Ptree).
        
        % Number of permutations (log) given the data structure.
        % It is assumed that Ptree was constructed using the
        % same X as input.
        lW = palm_maxshuf(Ptree,stype,true);
        varargout{1} = lW;
        
        % Number of permutations (log) if there were no data structure:
        if nargout > 1,
            lfac = palm_factorial(N);
            [~,~,S] = unique(X,'rows');
            U   = unique(S);
            nU  = numel(U);
            if strcmpi(stype,'perms') || strcmpi(stype,'both'),
                cnt = zeros(nU,1);
                for u = 1:nU,
                    cnt(u) = sum(S == U(u));
                end
                plW0 = lfac(N+1) - sum(lfac(cnt+1));
            else
                plW0 = 0;
            end
            if strcmpi(stype,'flips') || strcmpi(stype,'both'),
                cnt = zeros(nU,1);
                for u = 1:nU,
                    cnt(u) = sum(S == U(u));
                end
                slW0 = nU*log(2);
            else
                slW0 = 0;
            end
            lW0 = plW0 + slW0;
            varargout{2} = lW0;
        end
        
        % If the user wants, output also the Huberman & Hogg complexity,
        % which is computed recursively below
        if nargout > 2,
            varargout{3} = hhcomplexity(Ptree,1) - 1;
        end
        
    case 'distances',
        
        % Average change per permutation, i.e., average
        % Hamming distance.
        varargout{1} = mean(sum(bsxfun(@ne,Pset(:,1),Pset),1),2);
        
        % Average Euclidean distance per permutation.
        varargout{3} = mean(sum(bsxfun(@minus,Pset(:,1),Pset).^2,1).^.5,2);
        
        % For the Hamming and Euclidean, now take ties in X into account.
        % Also, compute the Spearman for each case
        if strcmpi(stype,'perms'),
            XP = X(Pset);
            varargout{5} = mean(1-6*sum(bsxfun(@minus,Pset(:,1),Pset).^2,1)/N/(N^2-1),2);
        elseif strcmpi(stype,'flips'),
            XP = bsxfun(@times,X,Pset);
            [~,iXP] = sort(XP);
            varargout{5} = mean(1-6*sum(bsxfun(@minus,iXP(:,1),iXP).^2,1)/N/(N^2-1),2);
        elseif strcmpi(stype,'both'),
            XP = X(abs(Pset));
            XP = sign(Pset).*XP;
            [~,iXP] = sort(XP);
            varargout{5} = mean(1-6*sum(bsxfun(@minus,iXP(:,1),iXP).^2,1)/N/(N^2-1),2);
        end
        varargout{2} = mean(sum(bsxfun(@ne,XP(:,1),XP),1),2);
        varargout{4} = mean(sum(bsxfun(@minus,XP(:,1),XP).^2,1).^.5,2);
end

% ==============================================================
function D = hhcomplexity(Ptree,D)
% Computes recursively the Huberman & Hogg complexity.
% For the 1st iteration, D = 1.

for u = 1:size(Ptree,1),
    if isnan(Ptree{u,1}(1)),
        k = size(Ptree{u,3},1);
    else
        k = numel(unique(Ptree{u,1}(:,1)));
    end
    D = D * (2^k - 1);
    if size(Ptree{u,3},2) > 1,
        D = hhcomplexity(Ptree{u,3},D);
    end
end
