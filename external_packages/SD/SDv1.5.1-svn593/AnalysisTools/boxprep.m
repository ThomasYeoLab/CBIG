% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================
function z = boxprep(varargin)% z = boxprep(varargin)% Prepare multiple inputs for boxplot call% A typical use is in conjunction with% the M-file boxplots. For example,% x = rand(20,4);% y = rand(1,30); % five sets of datda% names = strvcat('c1','c2','c3','c4','c5');% boxplots(boxprep(x,y),names,'random data')% June 16, November 12, 2003, Herman Gollwitzern = length(varargin);tmp = varargin{1};s = size(tmp);if s(1) == 1	z = tmp';else	z = tmp;endfor k = 2:n	tmp = varargin{k};	if size(tmp,1) == 1		tmp = tmp';	end	z = splicenan(z,tmp,2);endfunction z = splicenan(A,B,dim)%splice(A,B,dim) concatenates A and B even if the%sizes don't match up. Padding with nans%is carried out if data is numeric and%blanks if inputs are char.%Herman Gollwitzer, February 16, November 27, 1998% June 16, 2003. Use nan instead of zero for numeric arraysif nargin < 3	dim = 1;endif ~isequal(class(A),class(B))	error('Arrays must be of same class');end[m,n] = size(A);[p,q] = size(B);if ischar(A)	fill = blanks(1);else	fill = nan;endif dim == 1	z = repmat(fill,m+p,max(n,q));	z(1:m,1:n) = A; 	z(m+(1:p),1:q) = B;else	z = repmat(fill,max(m,p),n+q);	z(1:m,1:n) = A;	z(1:p,n+(1:q)) = B;end