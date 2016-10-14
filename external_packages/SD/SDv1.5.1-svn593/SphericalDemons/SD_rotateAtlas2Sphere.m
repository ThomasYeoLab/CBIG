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
function [curr_sbjWarp, Center1, Center2, Center3, prev_energy, curr_energy] = SD_rotateAtlas2Sphere(sbjMesh, basic_atlas, parms, SearchWidth, numIntervals)

% [curr_sbjWarp, Center1, Center2, Center3] = SD_rotateAtlas2Sphere(sbjMesh, basic_atlas, parms, SearchWidth)
% 
% Assume only 1 data dim!!

Center1 = 0; bestCenter1 = 0;
Center2 = 0; bestCenter2 = 0;
Center3 = 0; bestCenter3 = 0;

if(nargin < 5)
    numIntervals = 8;
end

if(nargin < 4)
    SearchWidth = 64/180*pi;
end

tmp_sbjWarp = parms.sbjWarp;

if(parms.rotate_latlon)
    if(parms.fast_rotation)
        num_per_side = floor(sqrt(size(basic_atlas.vertices, 2)/163842)*512);
        basic_atlas.im = MARS_convertMesh2Image(basic_atlas, [basic_atlas.mean; basic_atlas.variance], num_per_side, num_per_side);
    else
        basic_atlas.im = MARS_convertMesh2Image(basic_atlas, [basic_atlas.mean; basic_atlas.variance], 512, 512);
    end
end
curr_sbjWarp = tmp_sbjWarp;

curr_energy = inf;
while(SearchWidth >= pi/180)
    interval = 2*SearchWidth/numIntervals;
    for alpha = Center1-SearchWidth:interval:Center1+SearchWidth
        for beta = Center2-SearchWidth:interval:Center2+SearchWidth
            for gamma = Center3-SearchWidth:interval:Center3+SearchWidth

                tmp_sbjWarp.curr_vertices = MARS_zrotate(MARS_yrotate(MARS_zrotate(parms.sbjWarp.curr_vertices, gamma), beta), alpha);

                if(parms.rotate_latlon)
                    gauss_parms = MARS_bilinearInterpolate(tmp_sbjWarp.curr_vertices, basic_atlas.im);
                else
                    gauss_parms = MARS_linearInterpolate(tmp_sbjWarp.curr_vertices, basic_atlas, [basic_atlas.mean; basic_atlas.variance]);
                end
                
                tmp_energy = ComputeObjFnGivenData(sbjMesh.current_data, gauss_parms(1:end/2,:), gauss_parms(end/2+1:end,:));
                
                if(alpha == 0 && beta == 0 && gamma == 0)
                   prev_energy = tmp_energy; 
                end
                
                if(tmp_energy < curr_energy)
                    if(~parms.minimal_verbosity)
                        disp(['Rotate by ' num2str(alpha) ', ' num2str(beta) ', ' num2str(gamma)]);
                    end
                    curr_energy = tmp_energy;
                    curr_sbjWarp.curr_vertices = tmp_sbjWarp.curr_vertices;
                    bestCenter1 = alpha;
                    bestCenter2 = beta;
                    bestCenter3 = gamma;
                end


            end
        end
    end
    SearchWidth = SearchWidth/2;
    Center1 = bestCenter1;
    Center2 = bestCenter2;
    Center3 = bestCenter3;
end