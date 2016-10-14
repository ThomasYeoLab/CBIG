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
function new_vec = ParallelTransport(old_vec, old_pt, new_pt, radius)

% new_vector = ParallelTransport(old_vector, old_pt, new_pt, radius)
%
% Assumes everything is 3 x N

if(nargin < 4)
    radius = 100;
end

if(nargin < 5)
    tol = 1e-3;
end


if(max(abs(dot(old_vec, old_pt, 1)./radius./(sqrt(sum(old_vec.^2, 1))+eps))) > tol)
    keyboard;
    error(['ParallelTransport: old_direction should be perpendicular to old_pt: ' num2str(max(abs(dot(old_vec, old_pt, 1)./radius./sqrt(sum(old_vec.^2, 1)))))]); 
end

% note that oldw = neww
oldw = cross(old_pt, new_pt, 1);
oldw = oldw ./ repmat( sqrt(sum(oldw.^2, 1))+eps, 3, 1);

oldv = cross(oldw, old_pt, 1);
oldv = oldv ./ repmat( sqrt(sum(oldv.^2, 1)) +eps, 3, 1);

newv = cross(oldw, new_pt, 1);
newv = newv ./ repmat( sqrt(sum(newv.^2, 1)) +eps, 3, 1);

proj_on_newv = dot(old_vec, oldv, 1);
proj_on_neww = dot(old_vec, oldw, 1);

new_vec = repmat(proj_on_newv, 3, 1) .* newv + repmat(proj_on_neww, 3, 1) .* oldw;

% cleaning up for close points...
index = (sqrt(sum((old_pt - new_pt).^2, 1)) < 1e-4);
new_vec(:, index) = old_vec(:, index);

% clean up antipodal points (for some reason setting, to old vec for
% antipodal points rather than 0 leads to better stability.
index = (sqrt(sum((old_pt + new_pt).^2, 1)) < 1e-4);
new_vec(:, index) = old_vec(:, index);


% % This old implementation is outdated
% %
% % First compute great circle direction to new pt;
% oldv = MARS_findTangentVecPt1toPt2(old_pt, new_pt, radius);
% % if(abs(sum(oldv)) == 0)
% %     new_vec = old_vec;
% %     return;
% % end
% 
% oldv = oldv ./ repmat( sqrt(sum(oldv.^2, 1)) +eps, 3, 1);
% 
% 
% 
% oldw = cross(old_pt, oldv, 1);
% oldw = oldw ./ repmat( sqrt(sum(oldw.^2, 1))+eps, 3, 1);
% 
% % note that oldw = neww
% newv = cross(oldw, new_pt, 1);
% newv = newv ./ repmat( sqrt(sum(newv.^2, 1))+eps, 3, 1);
% 
% proj_on_newv = dot(old_vec, oldv, 1);  
% proj_on_neww = dot(old_vec, oldw, 1);
% 
% new_vec = repmat(proj_on_newv, 3, 1) .* newv + repmat(proj_on_neww, 3, 1) .* oldw;
% 
% % find points which were the same
% index = find(sum(abs(old_pt - new_pt), 1) == 0);
% new_vec(:, index) = old_vec(:, index);
% 
% % clean up antipodal points (for some reason setting, to old vec for
% % antipodal points rather than 0 leads to better stability.
% index = (sqrt(sum((old_pt + new_pt).^2, 1)) < 1e-4);
% new_vec(:, index) = old_vec(:, index);



