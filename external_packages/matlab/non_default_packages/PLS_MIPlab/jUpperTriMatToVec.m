function v=jUpperTriMatToVec(m,offset)
% converts the upper-triangular part of a matrix to a vector
%
% IN:
%   m: matrix
%   offset: offset above leading diagonal, fed to triu function
% OUT:
%   v: vector of upper-triangular values
% 
% v1.0 Oct 2009 Jonas Richiardi
% - initial release
% v2.0 Feb 2011 Jonas Richiardi
% - considerable speedup by vectorising,
% based on code by vortse a. Could also use Matlab FEX 23391
% by Bruno Luong


v=m(triu(true(size(m)),offset));