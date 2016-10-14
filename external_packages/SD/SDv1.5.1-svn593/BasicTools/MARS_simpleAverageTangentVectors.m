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
function smoothed_tangent_vectors = MARS_simpleAverageTangentVectors(uniform_mesh, tangent_vectors, var)

% smoothed_tangent_vectors = MARS_simpleAverageTangentVectors(uniform_mesh, tangent_vectors, var)
% 
% Written by Thomas, MIT
%
% Warning: We assume uniform_mesh to avoid for loop.
%
% Will not be using this method...

smoothed_tangent_vectors = tangent_vectors;

vertexNbors = uniform_mesh.vertexNbors;


if(size(vertexNbors, 1) ~= 6)
    error('Assumes uniform mesh with 6 neighbors!');
end


if(var > 0)
   factor = exp(-1/(2*var)); 
else
   factor = 1;
end

% neighbor 1
smoothed_tangent_vectors = smoothed_tangent_vectors + factor*ParallelTransport(tangent_vectors(:, vertexNbors(1, :)), uniform_mesh.vertices(:, vertexNbors(1, :)), uniform_mesh.vertices);

% neighbor 2
smoothed_tangent_vectors = smoothed_tangent_vectors + factor*ParallelTransport(tangent_vectors(:, vertexNbors(2, :)), uniform_mesh.vertices(:, vertexNbors(2, :)), uniform_mesh.vertices);

% neighbor 3
smoothed_tangent_vectors = smoothed_tangent_vectors + factor*ParallelTransport(tangent_vectors(:, vertexNbors(3, :)), uniform_mesh.vertices(:, vertexNbors(3, :)), uniform_mesh.vertices);

% neighbor 4
smoothed_tangent_vectors = smoothed_tangent_vectors + factor*ParallelTransport(tangent_vectors(:, vertexNbors(4, :)), uniform_mesh.vertices(:, vertexNbors(4, :)), uniform_mesh.vertices);

% neighbor 5
smoothed_tangent_vectors = smoothed_tangent_vectors + factor*ParallelTransport(tangent_vectors(:, vertexNbors(5, :)), uniform_mesh.vertices(:, vertexNbors(5, :)), uniform_mesh.vertices);

% neighbor 6
smoothed_tangent_vectors(:, vertexNbors(6, :)~=0) =  smoothed_tangent_vectors(:, vertexNbors(6, :)~=0) + ...
                factor*ParallelTransport(tangent_vectors(:, vertexNbors(6, vertexNbors(6, :)~=0)), uniform_mesh.vertices(:, vertexNbors(6, vertexNbors(6, :)~=0)), uniform_mesh.vertices(:, vertexNbors(6, :)~=0));

% Average
smoothed_tangent_vectors = smoothed_tangent_vectors ./ (repmat(sum(vertexNbors~=0, 1), 3, 1)*factor + 1);

