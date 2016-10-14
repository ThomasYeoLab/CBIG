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
function statMat = ComputeStatisticsOnDistMat(distMat, type)

% statMat = ComputeStatisticsOnDistMat(distMat, type)
% one type = 90 x 8
% if type == 1, then it's mean else, median

if(nargin < 2)
   type = 1; 
end

%load(distMatFilename);


num_structures = size(distMat, 2);
num_comparisons = size(distMat, 1);
num_subjects = (1 + sqrt(1 + 4*num_comparisons))/2; %compute number of subjects
disp(['Num_subjects: ' num2str(num_subjects)]);


statMat = zeros(num_structures, num_subjects);
z = zeros(1, num_subjects);
for i = 1:num_structures %for each structure
   

    struct_distMat = reshape(distMat(:, i), num_subjects-1, num_subjects);    
    
    if(type == 1)
    struct_distMat(struct_distMat == -1) = 0; % set missing structure comparison to 0.
    
    for j = 1:num_subjects
       z(j) = length(find(struct_distMat(:, j) ~= 0));
       
       if(z(j) == 0)
          z(j) = num_subjects - 1; %This is to prevent division by 0.
       end
    end
    
    statMat(i, :) = sum(struct_distMat, 1)./z;
    
    else
        
        statMat(i, :) = MedianIgnoringZeros(struct_distMat, 1)';
    end
    
end