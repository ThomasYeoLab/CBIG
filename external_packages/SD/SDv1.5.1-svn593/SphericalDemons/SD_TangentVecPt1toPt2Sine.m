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
function tangent_vec = SD_TangentVecPt1toPt2Sine(Pt1, Pt2, radius)

% Find tangent vectors from Pt1 to Pt2. Length determines by sine of angle
% between Pt1 to Pt2
%Pt1 is 3 x N
%Pt2 is 3 x N

if(nargin < 3)
    radius = 100;
end
    
debug = 1;

if(debug)
    total_length1 = sqrt(sum(Pt1.^2, 1));
    total_length2 = sqrt(sum(Pt2.^2, 1));
    
    if(max(abs(total_length1 - total_length2)) > 1e-4)
        disp(max(abs(total_length1 - total_length2)))
        warning('MARS_findTangentVecPt1toPt2: Pt1 and Pt2 not of same length');
    end
   
    Pt1 = Pt1./repmat(total_length1, 3, 1) * radius;
    Pt2 = Pt2./repmat(total_length2, 3, 1) * radius;
    
    total_length1 = sqrt(sum(Pt1.^2, 1));
    total_length2 = sqrt(sum(Pt2.^2, 1));
    
    if(max(abs(total_length1 - total_length2)) > 1e-4)
        disp(max(abs(total_length1 - total_length2)))
        error('MARS_findTangentVecPt1toPt2: Pt1 and Pt2 not of same length');
    end
        
end


%radius = sqrt(sum(Pt1(:,1).^2));

norm_vec = cross(Pt1, Pt2); %find normal vectors
tangent_vec = cross(norm_vec, Pt1); %find tangent Vectors

tangent_vec = tangent_vec /radius/radius; % normalize by radius^2 


%length_tangent_vec = sqrt(sum(tangent_vec.^2, 1));

%half_distance_between_pts =  sqrt(sum((Pt1-Pt2).^2, 1))/2;
%half_angle = asin(half_distance_between_pts./radius);

%correct_length = 2*half_angle*radius;

%index = find(length_tangent_vec ~= 0);

%if(~isempty(index))
%    tangent_vec(:,index) = tangent_vec(:,index)./repmat(length_tangent_vec(index), 3, 1) .* repmat(correct_length(index), 3,1);
%end

if(max(abs(dot(tangent_vec, Pt1, 1)./radius./(sqrt(sum(tangent_vec.^2, 1)) + eps))) > 1e-5)
    keyboard;
    error(['ParallelTransport: old_direction should be perpendicular to old_pt: ' num2str(max(abs(dot(tangent_vec, Pt1, 1)./radius./sqrt(sum(tangent_vec.^2, 1)))))]); 
end








