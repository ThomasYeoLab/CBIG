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
function WarpedPt = SD_warpPointbyTangentVecSine(Pt, tangentVec, input_radius)

%Assume Pt is 3 x N, tangentVec is 3 x N
%assume each of the point has the same radius.
%length of tangentVec corresponds to arclength

if(nargin <3)
    input_radius = 100;
end

%calculate the geodesic distance that has to be travelled by the warp
sine_angle = sqrt(sum(tangentVec.^2, 1))./input_radius;

%calculate distance that has to be traveled on the tangent plane at the
%starting point, before projecting onto the sphere to get destination
%points. (radius * tan \theta), but since we have sine(theta)
distance = input_radius .* sine_angle./sqrt(1 - sine_angle.^2);

%Find the set of tangentVec that are not equal to 0.
index = find(distance~=0);
WarpedPt = Pt;

% Move along tangent plane along tangent vector for a length specified by
% the variable distance
WarpedPt(:, index) = Pt(:, index) + tangentVec(:, index)./(repmat(sine_angle(:,index), 3, 1)*input_radius ).* repmat(distance(:,index), 3, 1);

% Now project onto the sphere of specified radius
new_radius = sqrt(sum(WarpedPt.^2, 1));
WarpedPt = WarpedPt./repmat(new_radius, 3, 1) * input_radius;




