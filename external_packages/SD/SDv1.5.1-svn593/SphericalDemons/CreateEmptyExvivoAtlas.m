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
function SD_atlas = CreateEmptyExvivoAtlas(hemi)

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

[status, hostname] = system('hostname -d');

if(~isempty(strfind(hostname, 'csail')))
    SD_atlas.parms.SUBJECTS_DIR = '/afs/csail.mit.edu/group/vision/medical-restricted/3/ythomas/ex_vivo_hist/';
    SD_atlas.parms.uniform_mesh_dir = '/afs/csail/group/vision/cortex/ythomas/work/MARS/';
else
    SD_atlas.parms.SUBJECTS_DIR = '/autofs/space/tensor_020/users/ythomas/ex_vivo_hist/';
    SD_atlas.parms.uniform_mesh_dir = '/autofs/space/tensor_017/users/thomas/work/MARS/';
end
%SD_atlas.parms.warp_filename = '';
%SD_atlas.parms.annot_filename = '';

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
SD_atlas.parms.flipFacesBool = 1;
SD_atlas.parms.WORK_DIR = 'SD';
SD_atlas.parms.uniform_meshes = {'ic4.tri', 'ic5.tri', 'ic6.tri', 'ic7.tri'};
SD_atlas.parms.subject_cell = {'pm14686'    'pm1696'    'pm18992'    'pm20784'    'pm28193'    'pm295'    'pm38281'    'pm54491'    'pm5694'    'pm6895'};
SD_atlas.parms.structures_of_interest = {'Unknown'    'V1'    'V2'    'BA44'    'BA45'    'BA4a'    'BA4p'    'BA6'    'BA2' 'MT'};