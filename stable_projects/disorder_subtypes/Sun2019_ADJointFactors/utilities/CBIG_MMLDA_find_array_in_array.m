function idx = CBIG_MMLDA_find_array_in_array(A, B)
% idx = CBIG_MMLDA_find_array_in_array(A, B)
%
% Find index of array A in array B in order. This function only works for vectors.
% 
% Input:
%   - A   : input array A with shorter length 
%   - B   : target array B with larger length
%
% Output:
%   - idx           : index array 
%
% Example:
%   idx = CBIG_MMLDA_find_array_in_array([1 2], [3 2 1]);
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

idx = arrayfun(@(x) find(B==x,1), A);