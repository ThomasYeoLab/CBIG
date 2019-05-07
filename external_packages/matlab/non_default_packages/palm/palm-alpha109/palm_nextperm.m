function [a,succ] = palm_nextperm(a)
% Given a sequence of integers "a", returns the next lexicographic
% permutation of this sequence. If "a" is already the last possible
% permutation, returns a vector of zeros of size(a).
% Note that to shuffle vectors, they must be supplied as
% column vectors (N by 1).
% 
% Usage:
% [a1,succ] = palm_nextperm(a)
% 
% a    : 2D array to be shuffled. Only the 1st column is
%        considered for the permutations. The rows as a whole
%        are shuffled together.
% a1   : Permuted sequence of values that corresponds to the
%        next lexicographic permutation.
% succ : If a is already the last possible permutation,
%        a1 = flipud(a) and succ is false.
%        Otherwise sucs is true.
% 
% This function is an implementation of the "Algorithm L",
% by D. Knuth (see "The Art of Computer Programming", Vol.4,
% Fasc.2: Generating All Tuples and Permutations.
% See also palm_algol.m to produce all possible permutations for
% a given sequence in a single function.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Feb/2012 (first version)
% Oct/2013 (this version)
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

% Algorithm L
% Step L2
n = size(a,1);
j = n - 1;
while j > 0 && a(j,1) >= a(j+1,1),
    j = j - 1;
end

% If this isn't yet the last permutation, bring up the next one.
if j > 0,
    
    % Step L3
    l = n;
    while a(j,1) >= a(l,1),
        l = l - 1;
    end
    tmp  = a(j,:);
    a(j,:) = a(l,:);
    a(l,:) = tmp;
    
    % Step L4
    k = j + 1;
    l = n;
    while k < l,
        tmp  = a(k,:);
        a(k,:) = a(l,:);
        a(l,:) = tmp;
        k = k + 1;
        l = l - 1;
    end
    
    % Was the permutation successful?
    succ = true;
    
else
    % If the input is the last permutation, then there is no next.
    % Return then the first shuffle and a successful flag "false"
    % that can be tested outside.
    a    = flipud(a);
    succ = false;
end