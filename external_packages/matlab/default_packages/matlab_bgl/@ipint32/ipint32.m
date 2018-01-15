function ipi = ipint32(a)
% IPDOUBLE Create a int32 type that supports inplace modification
% ipi = ipint32(v) creates an inplace int32 object from any int32
% matrix v.
%
% Example:
%    ipd = ipint32(ones(5));

ip = inplace(a);
ipi = class(struct([]),'ipint32',ip);

