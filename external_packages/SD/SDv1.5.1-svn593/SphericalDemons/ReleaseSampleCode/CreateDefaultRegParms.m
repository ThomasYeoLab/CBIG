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
function parms = CreateDefaultRegParms


% parms = CreateDefaultRegParms
%
% Written by Thomas Yeo, MIT
%
% parms.meshes = lists of multiresolution cell. Minimum should
% contain one level
% parms.noise is the noise variance estimate
% parms.max_step is the regularization that determines the step size to
% take in Newton's method
% parms.sbjWarp = previous sbjWarp

parms.radius = 100;

parms.meshes = cell(0); %e.g. is ic4,ic5,ic6,subject_mesh
parms.curr_level = 0; % This indexs the meshes
parms.iter = [];

parms.max_step = 2; %this is sigma_x
parms.sbjWarp = [];
parms.movie_flag = 0; % if movie_flag is on, will take snapshot.

parms.fast_deformation_smoothing = 1; % if 1 will do approximate speedup
parms.fast_deformation_level = 3; % switch on speed up at mesh level 3 and above

parms.smooth_displacement_var = [2 2 2 2]; %var for smoothing, followed by num_iters 
parms.smooth_displacement_iter = [10 10 10 10];
 
parms.smooth_velocity = 0;
parms.fast_velocity_smoothing = 1;
parms.fast_velocity_level = 3;
parms.smooth_velocity_var = [2 2 2 2];
parms.smooth_velocity_iter = [4 4 4 4];

parms.rotate = 1;
parms.rotate_latlon = 1; % if 1, use latlon grid for interpolation in rotation registration
parms.fast_rotation = 1; % if 1, then will do approximate speed up, by using less samples for smaller uniform meshes
parms.rotate_width = 64/180*pi;
parms.rotate_interval = 8;
parms.rotate_multiscale = 1;

parms.latlon_interpolation = 1; % if 1, then will project vertices on lat-lon grid. This will speed up the interpolation
parms.latlon_tight = 1; % if 1, then will use less samples for smaller uniform meshes
parms.latlon_level = 3; % level at which we switch to latlon interpolation.


parms.rkhs = 1;
parms.geodesic = 0; %Note that geodesic is not rigorous since dist(s,c) would be incompatible in the first and second optimization steps of demons
parms.approxH = 0;

parms.use_prev_warp = 0;
parms.final_unfold = 0;

% Possible options:
% verbose = 1, minimal_verbosity = 0 => print a lot
% verbose = 0, minimal_verbosity = 0 => print less
% verbose = 0, minimal_verbosity = 1 => print the least
% verbose = 1, minimal_verbosity = 1 => does not make sense
parms.verbose = 0;               
parms.minimal_verbosity = 1; 


