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
function vals = MARS_bilinearInterpolateDynamicImage(points, uniform_mesh, data)

% vals = MARS_bilinearInterpolateDynamicImage(points, uniform_mesh, data)
%
% Written by Thomas Yeo
%
% Here the input data is quickly projected onto latlon grid. Assumption is
% that uniform_mesh.numRows, uniform_mesh.numCols, uniform_mesh.weights,
% uniform_mesh.vno exists.
% 
% note that data "exists" on the vertices of uniform_mesh

num_data = size(data, 1);
im = repmat(uniform_mesh.weights(1, :), num_data, 1) .* data(:, uniform_mesh.vno(1,:)) + ...
         repmat(uniform_mesh.weights(2, :), num_data, 1) .* data(:, uniform_mesh.vno(2,:)) + ...
         repmat(uniform_mesh.weights(3, :), num_data, 1) .* data(:, uniform_mesh.vno(3,:));

im = reshape(im, num_data, uniform_mesh.numRows, uniform_mesh.numCols);
     
vals = MARS_bilinearInterpolate(points, im);




