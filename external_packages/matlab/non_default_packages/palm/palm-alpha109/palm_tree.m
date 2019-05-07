function Ptree = palm_tree(B,M)
% Generates a tree that takes the dependence structure
% between observations into account. From this tree, the
% permutations can be generated later.
%
% Usage:
% Ptree = palm_tree(B,M)
%
% - B       : Multi-level block definitions.
% - M       : Design matrix.
% - Ptree   : Permutation tree, from which permutations are generated
%             later.
%
% Each node is a cell with 4 elements:
% N{1,1}  : A 3-column array for whole block in which:
%           - the 1st is a sequence of indices that indicates the
%             current lexicographic permutation.
%           - the 2nd are indices that indicate the current
%             shuffling in relation to the original
%           - the 3rd are indices that indicate the current
%             permutation in relation to the previous
%          For within-block, this is a NaN.
% N{1,2} : A logical vector indicating the current state of sign
%          flips for the tree. 0 is treated as 1, and 1 as -1.
% N{1,3} : The branches that begin here.
%
% Reference:
% * Winkler AM, Webster MA, Vidaurre D, Nichols TE, Smith SM.
%   Multi-level block permutation. Neuroimage. 2015;123:253-68.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2013
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

if nargin == 1 || isempty(M),
    M = (1:size(B,1))';
elseif size(B,1) ~= size(M,1),
    error('The two inputs must have the same number of rows.');
end

% Make some initial sanity checks:
Bs = sortrows(B);
warned = checkblk(Bs,Bs(1)>=0,0);
if warned,
    error([...
        'Due to one or more of the issues listed above, the block\n' ...
        'definition may cause problems when the permutations are generated.\n'  ...
        'Please, correct the block definition file and try again.%s'],'');
end

% Order of the observations in the original data. Note
% that B should have not been sorted, and should retain
% the order as originally entered by the user (even if
% the rows are totally scrambled).
O = (1:size(M,1))';

% Now make the actual tree.
% The sanity of the block definitions should have already
% been taken care of by the wrapper function, so no need to
% check validity here. If the first element is negative,
% it is fair to assume that so are the remaining of the
% elements in the first column.
wholeblock = B(1) > 0;
Ptree = cell(1,3);
[Ptree{1},Ptree{3}] = maketree( ...
    B(:,2:end),M,O,wholeblock,wholeblock);
if wholeblock,
    Ptree{2} = false(size(Ptree{3},1),1);
else
    Ptree{2} = [];
end

% ==============================================================
function [S,Ptree] = maketree(B,M,O,wholeblock,nosf)
% Now makes the actual tree, each branch recursively.
% - B    : Block definitions
% - M    : Design matrix
% - O    : Observation indices
% - wholeblock : boolean, indicates if this is part of whole-block
%          at the immediately upper level.
% - nosf : It also implies no signflip this level, because
%          of a higher whole-block)

% Unique blocks at this level & some other vars for later
B1 = B(:,1);
U  = unique(B1);
nU = numel(U);
if size(B,2) > 1,
    Ptree = cell(nU,3);
else
    Ptree = cell(nU,1);
end

% For each block
for u = 1:nU,
    
    % Enter into each unique block
    idx = B1 == U(u);
    
    % If this isn't the last level, continue constructing
    % the branches recursively.
    if size(B,2) > 1,
        wholeblockb = B(find(idx,1),1) > 0;
        [Ptree{u,1},Ptree{u,3}] = ...
            maketree(             ...
            B(idx,2:end),         ...
            M(idx,:),             ...
            O(idx),               ...
            wholeblockb,          ...
            wholeblockb || nosf);
        Ptree{u,2} = [];
        
        % Count the number of possible sign-flips for these branches
        if nosf,
            % If it was flipped at higher levels (whole-block)
            Ptree{u,2} = [];
            
        elseif size(Ptree{u,3},2) > 1,
            % If it might be flipped here, but there are distal branches:
            % If this is whole-block, assign a number. If within-block,
            % no sign-flips allowed at this level.
            if isnan(Ptree{u,1}(1)),
                Ptree{u,2} = [];
            else
                Ptree{u,2} = false(size(Ptree{u,3},1),1);
            end            
        else
            % If there are no further branches
            Ptree{u,2} = false(size(Ptree{u,3},1),1);
        end
        
    else
        % At the terminal branches, there is no more tree, so
        % just keep track of the observation indices.
        Ptree{u,1} = O(idx);
    end
end

% Make the sequence of values that are the reference for the
% lexicographic permutations to be done later. This is doesn't
% apply to the highest level.
if wholeblock && nU > 1,
    
    % Identify repeated sets of rows, which receive then
    % the same index; these repetitions are taken care of
    % later by the Algorithm L.
    B1M     = sortrows([B1 M]); % B1 is here to be together with M during sorting
    Ms      = B1M(:,2:end);     % but its removed here, as it's of no use.
    [~,~,S] = unique(reshape(Ms',[numel(Ms)/nU nU])','rows');
    
    % Put in ascending order, and (un)shuffle the branches
    % accordingly
    [S,idx] = sort(S);
    S = [S (1:numel(S))' (1:numel(S))'];
    Ptree = Ptree(idx,:);
    
elseif wholeblock && nU == 1,
    % For whole block starting at the second level, the
    % permutation matrix is simply the number 1.
    S = [1 1 1];
    
else
    % There isn't whole-block permutation at the 1st level,
    % only within-block, so this case is marked as NaN.
    S = NaN;
end

% ==============================================================
function warned = checkblk(B,wholeblock,recdepth)
% For a given block, check if:
% - the leftmost column is valid.
% - the blocks are of the same size for whole-block permutation.
% - the sign indicator is the same for all blocks for whole-block
%   permutation.
% - the indices aren't 0 or non-integer.
% - the tree branches are all identical for whole-block permutation.
% Note this uses recursion.

% Vars for later
warned = false;
B1     = B(:,1);
U      = unique(B1);
nU     = numel(U);
Ucnt   = zeros(nU,1);
Usgn   = Ucnt;
Uvec   = cell(nU,1);

% Using positive/negative indices implies that the leftmost column must
% exist and be filled entirely with the same digit.
if recdepth == 0 && numel(U) > 1,
    error('The highest level (leftmost column) must be entirely filled by the same value.');
end

% For each block
for u = 1:nU,
    
    if U(u) == 0,
        
        % The index 0 cannot be considered valid (neither
        % positive or negative).
        warning('The digit 0 should be avoided as block indicator (level %d).\n',recdepth); %#ok
        warned = true;
        
    elseif rem(U(u),1) ~= 0,
        
        % Let's avoid fractional indices too.
        warning('Non-integer indices should be avoided (level %d).\n',recdepth); %#ok
        warned = true;
        
    else
        
        % Enter into each unique block to see what else
        idx = B1 == U(u);
        if size(B,2) > 1,
            
            % Here the test for whole-block permutation allows 0 as index,
            % just so that these otherwise invalid blocks are not ignored.
            % The warning message above should raise the user attention.
            checkblk(B(idx,2:end),B(find(idx,1),1)>=0,recdepth+1);
            
            % This is to check if the remainder cols are all equal. It's
            % necessary that the rows of the original B are sorted for this
            % to work.
            tmp = B(idx,2:end);
            Uvec{u} = tmp(:)';
        end
        
        % Check the size and sign of the sub-blocks
        Ucnt(u) = sum(idx);
        Usgn(u) = sign(U(u));
    end
end

% Check if all subblocks within each EB are of the same size.
if wholeblock && any(diff(Ucnt)),
    
    % Note that counting the number of times Matlab/Octave
    % reports as the source of the error the line above where
    % palm_renumber is called again also helps to identify
    % which is the offending block level.
    error('Not all sub-blocks within an EB are of the same size at level %d.\n',recdepth);
    
elseif wholeblock && numel(unique(Usgn)) > 1,
    
    % Make sure that for whole-block permutation, there is no mix of
    % sub-blocks with positive and negative signs, as these cannot be
    % shuffled.
    error('Not all sub-blocks within an EB have the same sign indicator at level %d.\n',recdepth);
    
elseif wholeblock,
    
    % Now check whether the branches that begin at a given level are all
    % identical, as needed for whole-block permutation. Note that for this
    % check to work, the blocks must be numbered using the "restart" (i.e.,
    % not "continuous") method. See the function palm_reindex.m.
    Uvec = cell2mat(Uvec);
    Ur = unique(Uvec,'rows');
    if size(Ur,1) > 1,
        warning('Not all blocks are identical after level %d to allow valid whole-block permutation.\n',recdepth); %#ok
        warned = true;
    end
end

