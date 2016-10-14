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
function [dice_vec, dice_overall] = MARS_computeDice(MARS_sbjMesh, OutputLabels, weightBool, weights)

% dice_vec = MARS_computeDice(MARS_sbjMesh, OutputLabels, weightBool, weights)

if(nargin < 4)
   weightBool = 0;
end

if(weightBool == 0)
   weights = ones(size(OutputLabels));
end

numLabels = MARS_sbjMesh.MARS_ct.numEntries;

dice_vec = zeros(numLabels, 1);
OutputLabels = int32(OutputLabels);

for i = 1:numLabels

    gt_struct = MARS_sbjMesh.MARS_label == int32(i);
    out_struct = OutputLabels == int32(i);

    overlap = sum(double(gt_struct.*out_struct).*weights);
    gt_area = sum(double(gt_struct).*weights);
    out_area = sum(double(out_struct).*weights);

    dice_vec(i) = 2*overlap/(gt_area+out_area);
end


overlap = sum(double((OutputLabels == MARS_sbjMesh.MARS_label)).*weights);
area = sum(weights);

dice_overall = overlap/area;