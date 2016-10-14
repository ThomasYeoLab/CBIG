function ipd = ipdouble(a)
% IPDOUBLE Create a double type that supports inplace modification
% ipd = ipdouble(v) creates an inplace double object from any double
% matrix v.
%
% Example:
%    ipd = ipdouble(ones(5));

ip = inplace(a);
ipd = class(struct([]),'ipdouble',ip);

