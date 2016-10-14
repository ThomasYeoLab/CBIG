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
function distMat = MARS_computeDistBetAnnotationAllSbj(MARS_sbjMeshCell, structures, subject_index)

% distMat = MARS_computeDistBetAnnotationAllSbj(MARS_sbjMeshCell, structures, subject_index)
% if subject_index is given. Only compare SbjMeshCell{subject_index} with
% all the other subjects.

num_subjects = length(MARS_sbjMeshCell);

if(nargin < 3)
    distMat = zeros(num_subjects*(num_subjects-1), length(structures));
    count = 0;
    for i = 1:num_subjects
        for j = 1:num_subjects
            if(i~=j)
                disp([num2str(i) ',' num2str(j)]);
                tic
                count = count + 1;
                distMat(count, :) = MARS_computeMeanDistanceBetAnnotation(MARS_sbjMeshCell{i}, MARS_sbjMeshCell{j}, structures);
                toc
            end
        end
    end
else
    distMat = zeros(2*(num_subjects-1), length(structures));
    
    count = 0;
    j = subject_index;
    for i = 1:num_subjects
        if(i~=j)
                disp([num2str(i) ',' num2str(j)]);
                tic
                count = count + 1;
                distMat(count, :) = MARS_computeMeanDistanceBetAnnotation(MARS_sbjMeshCell{i}, MARS_sbjMeshCell{j}, structures);
                count = count + 1;
                distMat(count, :) = MARS_computeMeanDistanceBetAnnotation(MARS_sbjMeshCell{j}, MARS_sbjMeshCell{i}, structures);
                toc
        end  
    end
end