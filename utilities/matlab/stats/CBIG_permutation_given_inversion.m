function permut = CBIG_permutation_given_inversion(inv_vec)

% Given the inversion of a permutation (4 3 2 1 5), get the permutation.
%
%   permut = CBIG_permutation_given_inversion(inv_vec)
%   Input:
%       inv_vec : inversion vector e.g. 4 3 2 1 5
%   Output:
%       permut  : permutation e.g. 5 4 3 2 1
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

n = length(inv_vec);
permut = n;

for i = n-1:-1:1
    permut = [permut(1:inv_vec(i)) i permut(inv_vec(i)+1:end)];
end
