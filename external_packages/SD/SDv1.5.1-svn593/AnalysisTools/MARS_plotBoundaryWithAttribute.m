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
function s = MARS_plotBoundaryWithAttribute(mesh, data, labels, threshold2)

% MARS_plotBoundaryWithAttribute(mesh, data, labels, threshold2)

BoundaryVec = zeros(length(labels), 1);
maxNeighbors = size(mesh.vertexNbors, 1);

for i = 1:length(labels)
    
    label_vertex = int32(labels(i));
    
    for k = 1:maxNeighbors
        v_neighbor = mesh.vertexNbors(k, i);
        if(v_neighbor ~= 0 && int32(labels(v_neighbor)) ~= label_vertex)
            BoundaryVec(i) = 1;
        end
    end
end

if(nargin < 4)
    data(BoundaryVec == 1) = max(data);
else
    threshold1 = max(data) - threshold2/10;
    data((BoundaryVec == 1)' & data > threshold2) = threshold1;
end

s = TrisurfMeshData(mesh, data, 1);

if(nargin > 3)
    numColors = 2000;
    
    upper_range = ceil((max(data)-threshold1)/(max(data)-min(data))*numColors);
    lower_range = ceil((threshold2 - min(data))/(max(data)-min(data)) * numColors);    
    
    buffer = round(0.25*(upper_range+lower_range));
    z = jet(lower_range+upper_range+buffer);
    
    new_colormap = zeros(numColors, 3);
    new_colormap(1:lower_range, :) = z(1:lower_range,:); 
    new_colormap(end-upper_range+1:end, :) = z(lower_range+buffer +1:end,:); 
    
    
    colormap(new_colormap);
    
end

shading interp;