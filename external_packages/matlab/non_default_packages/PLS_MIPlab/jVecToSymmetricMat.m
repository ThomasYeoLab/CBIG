function m=jVecToSymmetricMat(v,msize,offset)
% converts a vector into a square symmetric matrix
%
% IN:
%   v:      nx1 vector of numbers
%   msize:  scalar, resulting matrix size (msize rows by msize columns)
% OUT:
%   m: square symmetric matrix
% 
% v1.0 Nov 2009 Jonas Richiardi
% - initial release
% v1.0.1 Nov 2012 JR
% - fixed doc
% v2 Jan 2018 VK
% - added offset input in case diagonal is not full of zeros

v=v(:);

% get indices of upper triangular part (Peter Acklam's trick)
[m_i m_j] = find(triu(ones(msize),offset));

if (numel(m_i) ~= numel(v))
    error('final matrix size does not correspond to number of element in vector');
end

% copy vector to upper-tri matrix
m=zeros(msize);
for v_idx=1:numel(m_i)
    m(m_i(v_idx),m_j(v_idx))=v(v_idx);
end
% make symmetric
for r=1:msize(1)
    for c=1:r
        m(r,c)=m(c,r);
    end
end
