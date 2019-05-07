function [unsrtR,S,srtR] = palm_competitive(X,ord,mod)
% Sort a set of values and return their competition
% ranks, i.e., 1224, or the modified competition ranks,
% i.e. 1334. This makes difference only when there are
% ties in the data. The function returns the ranks in
% their original order as well as sorted.
% 
% Usage:
% [unsrtR,S,srtR] = palm_competitive(X,ord,mod)
% 
% Inputs:
% - X      : 2D array with the original data. The
%            function operates on columns. To operate
%            on rows or other dimensions, use transpose
%            or permute the array's higher dimensions.
% - ord    : Sort as 'ascend' (default) or 'descend'.
% - mod    : If true, returns the modified competition
%            ranks, i.e., 1334. This is the
%            correct for p-values and cdf. Otherwise
%            returns standard competition ranks.
% 
% Outputs:
% - unsrtR : Competitive ranks in the original order.
% - S      : Sorted values, just as in 'sort'.
% - srtR   : Competitive ranks sorted as in S.
%
% Examples:
% - To obtain the empirical cdf of a dataset in X, use:
%   cdf   = palm_competitive(X,'ascend',true)/size(X,1);
% - To obtain the empirical p-values for each value in X, use:
%   pvals = palm_competitive(X,'descend',true)/size(X,1);
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Nov/2012
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

% Check inputs
if nargin < 1 || nargin > 3,
    error('Incorrect number of arguments.');
elseif nargin == 1,
    ord = 'ascend';
    mod = false;
elseif nargin == 2,
    mod = false;
end

% Important: The function starts by computing the
% unmodified competition ranking. This can be
% ascending or descending. However, note that the
% modified ranking for the ascending uses the
% unmodified descending, whereas the modified
% descending uses the modified ascending, hence
% the need to "reverse" the inputs below.
if mod,
    if strcmpi(ord,'ascend'),
        ord = 'descend';
    elseif strcmpi(ord,'descend'),
        ord = 'ascend';
    end
end

% Unmodified competition ranking
[nR,nC] = size(X);
unsrtR  = single(zeros(size(X)));
[S,tmp] = sort(X,ord);
[~,rev] = sort(tmp);
srtR = repmat((1:nR)',[1 nC]);
for c = 1:nC, % loop over columns

    % Check for +Inf and -Inf and replace them
    % for a value just higher or smaller than
    % the max or min, respectively.
    infpos = isinf(S(:,c)) & S(:,c) > 0;
    infneg = isinf(S(:,c)) & S(:,c) < 0;
    if any(infpos),
        S(infpos,c) = max(S(~infpos,c)) + 1;
    end
    if any(infneg),
        S(infneg,c) = min(S(~infneg,c)) - 1;
    end
    
    % Do the actual sorting, checking for obnoxious NaNs
    dd = diff(S(:,c));
    if any(isnan(dd)),
        error(['Data cannot be sorted. Check for NaNs that might be present,\n', ...
            'or precision issues that may cause over/underflow.\n', ...
            'If you are using "-approx tail", consider adding "-nouncorrected".%s'],'');
    end
    f = find([false; ~logical(dd)]);
    for pos = 1:numel(f),
        srtR(f(pos),c) = srtR(f(pos)-1,c);
    end
    unsrtR(:,c) = single(srtR(rev(:,c),c)); % original order as the data
    
    % Put the infinities back
    if any(infpos),
        S(infpos,c) = +Inf;
    end
    if any(infneg),
        S(infneg,c) = -Inf;
    end
end

% Prepare the outputs for the modified rankings, i.e.,
% flip the sorted values and ranks
if mod,
    
    % Do the actual modification
    unsrtR = nR - unsrtR + 1;
    
    % Flip outputs
    if nargout >= 2,
        S = flipud(S);
    end
    if nargout == 3,
        srtR = flipud(nR - srtR + 1);
    end
end
