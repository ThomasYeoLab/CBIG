function ipa = assign(ipa,y,varargin)
% INPLACE/ASSIGN Assign an entire inplace array.
%
% assign(ipa,y) overwrites the contents of ipa with y.  
%
% Example:
%    ipa = inplace(ones(5,1));
%    assign(ipa,pi*ones(5,1));

ipd.assign(y);