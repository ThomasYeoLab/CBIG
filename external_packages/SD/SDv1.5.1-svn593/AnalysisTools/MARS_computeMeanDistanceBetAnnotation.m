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
function distVec = MARS_computeMeanDistanceBetAnnotation(MARS_sbjMesh1, MARS_sbjMesh2, structures, sbjWarp1, sbjWarp2, MARS_uniformMesh_SpatialPrior)

% MARS_computeMeanDistanceBetAnnotation(MARS_sbjMesh1, MARS_sbjMesh2, structures)
% Projection is done onto MARS_sbjMesh1
%
% This includes an alternative version where MARS_sbjMesh1 and
% MARS_sbjMesh2 are not in the same space, in which case, they are related
% via the warps sbjWarp1 and sbjWarp2.

if(nargin <= 3)
    labels2 = MARS_NNInterpolate(MARS_sbjMesh1.vertices, MARS_sbjMesh2, MARS_sbjMesh2.MARS_label);
else
    MARS_uniformMesh_SpatialPrior.vertices = sbjWarp1.curr_vertices;
    
    labels2 = MARS_NNInterpolate(sbjWarp2.curr_vertices, MARS_sbjMesh2, MARS_sbjMesh2.MARS_label);
    labels2 = MARS_NNInterpolate(MARS_sbjMesh1.vertices, MARS_uniformMesh_SpatialPrior, labels2);
end

% unscaling...
MARS_sbjMesh1.metricVerts = MARS_sbjMesh1.metricVerts/MARS_sbjMesh1.surface_scaling_factor;
labels1 = MARS_sbjMesh1.MARS_label;


num_vertices = int32(size(MARS_sbjMesh1.vertices, 2));
maxNeighbors = int32(size(MARS_sbjMesh1.vertexNbors, 1));

vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh1.vertexNbors, 1)), int32(size(MARS_sbjMesh1.vertices, 2)), int32(MARS_sbjMesh1.vertexNbors), single(MARS_sbjMesh1.metricVerts));
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);


distVec = zeros(1, length(structures));
for i = 1:length(structures)

    k = structures(i);
    %disp(num2str(k));
    %tic
    if(~isempty(find(labels1==k)) && ~isempty(find(labels2 == k))) %structure exists for both subject
        
        objVec1 = (labels1 == k); 
        boundaryVec1 = ConvertObjVec2Boundary(MARS_sbjMesh1, objVec1);
        objVec2 = (labels2 == k);
        boundaryVec2 = ConvertObjVec2Boundary(MARS_sbjMesh1, objVec2);
        
        final_cost1 = MARS_distSrcs2Dests(int32(boundaryVec1), int32(boundaryVec2), int32(num_vertices), int32(maxNeighbors), int32(MARS_sbjMesh1.vertexNbors), double(vertexDist2Nbors));
        final_cost2 = MARS_distSrcs2Dests(int32(boundaryVec2), int32(boundaryVec1), int32(num_vertices), int32(maxNeighbors), int32(MARS_sbjMesh1.vertexNbors), double(vertexDist2Nbors));
        
        final_cost1 = final_cost1(logical(boundaryVec2));
        final_cost2 = final_cost2(logical(boundaryVec1));
        
        distVec(i) = mean([final_cost1; final_cost2]);
        
    else
        distVec(i) = -1;
    end
    %toc
end


