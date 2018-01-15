function d=edit_distance_weighted(s,t,weight_delete,weight_insert,weight_replace)
  % Calculate a weighted string distance. If the weights are 1,1,1, then
  % the calculated distance is equal to the Levenshtein distance.
  %
  % @author B. Schauerte <schauerte@ieee.org>
  % @date 2010
  
  % Copyright 2010 B. Schauerte. All rights reserved.
  % 
  % Redistribution and use in source and binary forms, with or without 
  % modification, are permitted provided that the following conditions are 
  % met:
  % 
  %    1. Redistributions of source code must retain the above copyright 
  %       notice, this list of conditions and the following disclaimer.
  % 
  %    2. Redistributions in binary form must reproduce the above copyright 
  %       notice, this list of conditions and the following disclaimer in 
  %       the documentation and/or other materials provided with the 
  %       distribution.
  % 
  % THIS SOFTWARE IS PROVIDED BY B. SCHAUERTE ''AS IS'' AND ANY EXPRESS OR 
  % IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
  % WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
  % DISCLAIMED. IN NO EVENT SHALL B. SCHAUERTE OR CONTRIBUTORS BE LIABLE 
  % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
  % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
  % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
  % BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
  % WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
  % OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
  % ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  % 
  % The views and conclusions contained in the software and documentation
  % are those of the authors and should not be interpreted as representing 
  % official policies, either expressed or implied, of B. Schauerte.
  
  m=numel(s);
  n=numel(t);
  
  if nargin < 3, weight_delete=1; end
  if nargin < 4, weight_insert=1; end
  if nargin < 5, weight_replace=1; end
  
  d=zeros(m+1,n+1);
  
  % initialize distance matrix
  for i=1:m % deletion
    d(i+1,1)=d(i,1) + weight_delete;
  end
  for j=1:n % insertion
    d(1,j+1)=d(1,j) + weight_insert;
  end
  
  for j=2:n+1
    for i=2:m+1
      if s(i-1) == t(j-1)
        cost_replace=0;
      else
        cost_replace=weight_replace;
      end
      d(i,j)=min([ ...
        d(i-1,j) + weight_insert, ... % insertion
        d(i,j-1) + weight_delete, ... % deletion
        d(i-1,j-1) + cost_replace ... % substitution
        ]);
    end
  end
  
  d=d(m+1,n+1);
