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
function MARS_plot3BoundaryWithAttribute(mesh, data, label1, label2, label3)

% MARS_plot2BoundaryWithAttribute(mesh, data, label1, label2)

BoundaryVec1 = zeros(length(label1), 1);
maxNeighbors = size(mesh.vertexNbors, 1);

for i = 1:length(label1)
    
    label_vertex = int32(label1(i));
    
    for k = 1:maxNeighbors
        v_neighbor = mesh.vertexNbors(k, i);
        if(v_neighbor ~= 0 && int32(label1(v_neighbor)) ~= label_vertex)
            BoundaryVec1(i) = 1;
        end
    end
end

BoundaryVec2 = zeros(length(label2), 1);
for i = 1:length(label2)
    
    label_vertex = int32(label2(i));
    
    for k = 1:maxNeighbors
        v_neighbor = mesh.vertexNbors(k, i);
        if(v_neighbor ~= 0 && int32(label2(v_neighbor)) ~= label_vertex)
            BoundaryVec2(i) = 1;
        end
    end
end

BoundaryVec3 = zeros(length(label3), 1);
for i = 1:length(label3)
    
    label_vertex = int32(label3(i));
    
    for k = 1:maxNeighbors
        v_neighbor = mesh.vertexNbors(k, i);
        if(v_neighbor ~= 0 && int32(label3(v_neighbor)) ~= label_vertex)
            BoundaryVec3(i) = 1;
        end
    end
end




data(BoundaryVec1 == 1) = max(data);
data(BoundaryVec2 == 1) = min(data);
data(BoundaryVec3 == 1) = (max(data)+min(data))/2;
TrisurfMeshData(mesh, data, 1);

