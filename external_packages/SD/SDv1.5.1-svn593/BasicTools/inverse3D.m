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
function mat3D = inverse3D(input_matrix)

% input is 3 x 3 x N
% output is a 3 x 3 x N matrix where entry of each 3x3 is the inverse

cofactors = cofactors_sym(input_matrix);

det3D = sum(squeeze(input_matrix(1, :, :)) .* squeeze(cofactors(1, :, :)), 1) + eps;
mat3D = shiftdim(repmat(1./det3D, [3 1 3]), 2) .* cofactors;


function mat3D = cofactors_sym(input_matrix)

% input a 3 x 3 x N
% output is a 3 x 3 x N matrix where entry of each 3x3 is the cofactors

mat3D = zeros(size(input_matrix));
mat3D(1, 1, :) = input_matrix(2, 2, :) .* input_matrix(3, 3, :) - input_matrix(2, 3, :) .* input_matrix(3, 2, :);
mat3D(1, 2, :) = -1 * (input_matrix(2, 1, :) .* input_matrix(3, 3, :) - input_matrix(3, 1, :) .* input_matrix(2, 3, :));
mat3D(1, 3, :) = input_matrix(2, 1, :) .* input_matrix(3, 2, :) - input_matrix(2, 2, :) .* input_matrix(3, 1, :);
mat3D(2, 1, :) = mat3D(1, 2, :);
mat3D(2, 2, :) = input_matrix(1, 1, :) .* input_matrix(3, 3, :) - input_matrix(1, 3, :) .* input_matrix(3, 1, :);
mat3D(2, 3, :) = -1 * (input_matrix(1, 1, :) .* input_matrix(3, 2, :) - input_matrix(1, 2, :) .* input_matrix(3, 1, :));
mat3D(3, 1, :) = mat3D(1, 3, :);
mat3D(3, 2, :) = mat3D(2, 3, :);
mat3D(3, 3, :) = input_matrix(1, 1, :) .* input_matrix(2, 2, :) - input_matrix(1, 2, :) .* input_matrix(2, 1, :);


