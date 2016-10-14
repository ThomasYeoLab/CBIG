function n=end(ipa,a,n)
% INPLACE/END Allow inline "end" indices to work.
%
% n=end(ipa,a,n) Returns the ending index along dimension a.
%
% Example:
%    ipa = inplace(cumsum(ones(5,1)));
%    ipa(3:end)

n = size(ipa.get_a(),a);