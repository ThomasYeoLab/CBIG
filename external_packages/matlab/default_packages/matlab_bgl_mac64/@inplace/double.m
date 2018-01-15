function d = double(ipa)
% INPLACE/DOUBLE Return the data as a Matlab double type.  
% d = double(ipa) returns a double array.  Changing d will not change the
% ipa data.
%
% Example:
%    ipa = inplace(cumsum(ones(5,1)));
%    d = double(ipa);
%    d(2) = 10: % Doesn't change ipa
%    ipa

d = ipa.get_a();