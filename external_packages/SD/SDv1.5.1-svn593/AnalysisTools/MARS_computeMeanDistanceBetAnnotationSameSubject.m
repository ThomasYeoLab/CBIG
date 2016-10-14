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
function [distVec, distVecOverall] = MARS_computeMeanDistanceBetAnnotationSameSubject(MARS_sbjMesh, labels1, labels2)

% distVec = MARS_computeMeanDistanceBetAnnotation(MARS_sbjMesh, labels1, labels2)

num_structures = MARS_sbjMesh.MARS_ct.numEntries;

vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)), int32(size(MARS_sbjMesh.vertices, 2)), int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.metricVerts));
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);

num_vertices = int32(size(MARS_sbjMesh.vertices, 2));
maxNeighbors = int32(size(MARS_sbjMesh.vertexNbors, 1));


distVec = zeros(num_structures, 1);
for k = 1:num_structures

    if(~isempty(find(labels1==k)) && ~isempty(find(labels2 == k))) %structure exists for both subject
        
        objVec1 = (labels1 == k); 
        boundaryVec1 = ConvertObjVec2Boundary(MARS_sbjMesh, objVec1);
        objVec2 = (labels2 == k);
        boundaryVec2 = ConvertObjVec2Boundary(MARS_sbjMesh, objVec2);
        
        final_cost1 = MARS_distSrcs2Dests(int32(boundaryVec1), int32(boundaryVec2), int32(num_vertices), int32(maxNeighbors), int32(MARS_sbjMesh.vertexNbors), double(vertexDist2Nbors));
        final_cost2 = MARS_distSrcs2Dests(int32(boundaryVec2), int32(boundaryVec1), int32(num_vertices), int32(maxNeighbors), int32(MARS_sbjMesh.vertexNbors), double(vertexDist2Nbors));
        
        final_cost1 = final_cost1(logical(boundaryVec2));
        final_cost2 = final_cost2(logical(boundaryVec1));
    
        distVec(k) = mean([final_cost1; final_cost2]);
        
    else
        distVec(k) = -1;
    end
    %toc
end

distVecOverall = mean(distVec(distVec~=-1));