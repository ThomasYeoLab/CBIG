function ip=subsasgn(ip,S,B)
% INPLACE/SUBSASGN Support Matlab subscript assignment
%
% This function supports subscript assignment just like Matlab matrices.
%
% Example:
%    ipa = inplace(cumsum(ones(5,1)));
%    ipa([2 4]) = 1;
%    ipa

ip.subsasgn(S,B);