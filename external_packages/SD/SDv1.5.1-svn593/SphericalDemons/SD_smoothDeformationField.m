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
function sbjWarp = SD_smoothDeformationField(parms, vec_bool)

% sbjWarp = MARS_smoothDeformationField(parms, vec_bool)
% 
% Written by Thomas Yeo, MIT

if(nargin < 2)
   vec_bool = 0; 
end


sbjWarp = parms.sbjWarp;    
distSq2Nbr = parms.meshes{parms.curr_level}.vertexDistSq2Nbors(1);

if(~vec_bool)
    if(parms.fast_deformation_smoothing && parms.curr_level >= parms.fast_deformation_level)
        
        if(parms.verbose)
            disp('Approximate Smoothing');
        end
        % Compute input radius
        input_radius = sqrt(sum(parms.meshes{parms.curr_level}.vertices(:,1).^2));

        % Smooth
        sbjWarp.curr_vertices = MARS_AverageData(parms.meshes{parms.curr_level}, sbjWarp.curr_vertices, parms.smooth_displacement_var(parms.curr_level)*distSq2Nbr, ...
            parms.smooth_displacement_iter(parms.curr_level));

        % Reproject as points on the sphere
        new_radius = sqrt(sum(sbjWarp.curr_vertices.^2, 1));
        sbjWarp.curr_vertices = sbjWarp.curr_vertices./repmat(new_radius, 3, 1) * input_radius;
    else
        if(parms.verbose)
            disp('More exact Smoothing');
        end
        sbjWarp.curr_vertices = MARS_AverageDisplacementVectors(parms.meshes{parms.curr_level}, sbjWarp.curr_vertices, ...
            parms.smooth_displacement_var(parms.curr_level), parms.smooth_displacement_iter(parms.curr_level));
    end
else
    if(parms.fast_velocity_smoothing && parms.curr_level >= parms.fast_velocity_level)
        
        if(parms.verbose)
            disp('Approximate Tangent Vector Smoothing');
        end
        for i = 1:parms.smooth_velocity_iter(parms.curr_level)
            sbjWarp.vec = MARS_simpleAverageData(parms.meshes{parms.curr_level}, sbjWarp.vec, parms.smooth_displacement_var(parms.curr_level)*distSq2Nbr);
            sbjWarp.vec = MARS_projectGradOntoTangentPlane(parms.meshes{parms.curr_level}.vertices, sbjWarp.vec);
        end
    else
        if(parms.verbose)
            disp('Exact Tangent Vector Smoothing');
        end
        for i = 1:parms.smooth_velocity_iter(parms.curr_level)
            sbjWarp.vec = MARS_simpleAverageTangentVectors(parms.meshes{parms.curr_level}, sbjWarp.vec, parms.smooth_velocity_var(parms.curr_level));
            sbjWarp.vec = MARS_projectGradOntoTangentPlane(parms.meshes{parms.curr_level}.vertices, sbjWarp.vec);
        end
    end
end
