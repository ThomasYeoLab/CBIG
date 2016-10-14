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
function SD_atlas = CreateDefaultAtlasParm(hemi)

SD_atlas.zrotate = pi/90;
SD_atlas.yrotate = pi/90;
SD_atlas.radius = 100;
SD_atlas.basic_atlas = [];

% These parameters are mostly used for reading subjects and meshes

% parameters to be filled
if(nargin == 1)
    SD_atlas.parms.hemi = hemi;
else
    SD_atlas.parms.hemi = 'lh';
end
SD_atlas.multidim = 0;
SD_atlas.parms.scale = 0;
SD_atlas.parms.speed = 0;

SD_atlas.parms.SUBJECTS_DIR = fullfile('..','..','example_surfaces');
SD_atlas.parms.uniform_mesh_dir = fullfile('..','..');

% parameters likely to stay the same
SD_atlas.relative_variance = 1; % 1 implies using relative variance for thresholding else uses absolute value.
SD_atlas.relative_variance_val = 0.0001; %if relative variance is 0, will use this absolute value for thresholding!
SD_atlas.accumulate_variance = 1; % if 1, then accumulates variance only lower resolution mesh and then interpolates onto higher resolution mesh.

SD_atlas.parms.read_surface = @MARS2_readSbjMesh;
SD_atlas.parms.radius = 100;
SD_atlas.parms.surf_filename = 'sphere';
SD_atlas.parms.data_filename_cell = {'inflated.H','sulc', 'sulc', 'curv'}; 
SD_atlas.parms.normalizeBool = 1;
SD_atlas.parms.unfoldBool = 1;
SD_atlas.parms.flipFacesBool = 1; % Probably do not want to change this!!
SD_atlas.parms.WORK_DIR = 'SD';
SD_atlas.parms.uniform_meshes = {'ic4.tri', 'ic5.tri', 'ic6.tri', 'ic7.tri'};
SD_atlas.parms.subject_cell = {'OAS1_0001_MR1' ,...
    'OAS1_0002_MR1' ,...
    'OAS1_0003_MR1' ,...
    'OAS1_0004_MR1' ,...
    'OAS1_0005_MR1'};
