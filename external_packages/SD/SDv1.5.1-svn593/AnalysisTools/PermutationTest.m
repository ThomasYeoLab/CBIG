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
function p_value = PermutationTest(vec1, vec2, N, type)

% p_value = PermutationTest(vec1, vec2)
% 
% Test for mean(vec1) > mean(vec2)
% vec1 not necessarily same length as vec2
% 
% N is number of permutations.

if(type)
    true_diff = mean(vec1) - mean(vec2);
else
    true_diff = median(vec1) - median(vec2); 
end
%disp(num2str(true_diff_mean));
rand('state',sum(100*clock))

Length1 = length(vec1);
Length2 = length(vec2);

if(size(vec1,1) > size(vec1,2))
   combinedVec = [vec1;vec2];
else
   combinedVec = [vec1 vec2]; 
end
TotalLength = length(combinedVec);


count = 0;
for i = 1:(N-1)
    
    newVec = combinedVec(randperm(TotalLength));
    
    if(type)
        newDiff = mean(newVec(1:Length1)) - mean(newVec(Length1+1:end));
    else
        newDiff = median(newVec(1:Length1)) - median(newVec(Length1+1:end)); 
    end
        
    if(newDiff > true_diff)
        count = count + 1;
    end
end


p_value = (count+1)/N;

