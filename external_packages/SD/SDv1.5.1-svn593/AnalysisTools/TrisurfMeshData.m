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
function s = TrisurfMeshData(mesh, data, new_figure)

% TrisurfMeshData(mesh, data)

if(nargin > 2)
    figure;set(gcf, 'renderer', 'zbuffer');
end
s = trisurf(mesh.faces', mesh.vertices(1,:)', mesh.vertices(2,:)', mesh.vertices(3,:)', data);
axis equal;
% shading flat;

% 'unknown '
%     'bankssts '
%     'caudalanteriorcingulate '
%     'caudalmiddlefrontal '
%     'corpuscallosum '
%     'cuneus '
%     'entorhinal '
%     'fusiform '
%     'inferiorparietal '
%     'inferiortemporal '
%     'isthmuscingulate '
%     'lateraloccipital '
%     'lateralorbitofrontal '
%     'lingual '
%     'medialorbitofrontal '
%     'middletemporal '
%     'parahippocampal '
%     'paracentral '
%     'parsopercularis '
%     'parsorbitalis '
%     'parstriangularis '
%     'pericalcarine '
%     'postcentral '
%     'posteriorcingulate '
%     'precentral '
%     'precuneus '
%     'rostralanteriorcingulate '
%     'rostralmiddlefrontal '
%     'superiorfrontal '
%     'superiorparietal '
%     'superiortemporal '
%     'supramarginal '
%     'frontalpole '
%     'temporalpole '
%     'transversetemporal '