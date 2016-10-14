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
% 
%
% Just type "CoregisterSurfaces" without any arguments and things should work!
% This function co-registers surfaces. There are 4 input parameters  
% and 6 set of parameters you can play with (default should work ok). More explanations below.  
%
% What to expect: 
%       By default, there are 5 subjects. Each subject registration takes about 5 mins, 
%       so one round of co-registration for ALL subjects should take approximately 30 min. 
%       Also by default, there will be 2 rounds of co-registration, so total time is approximately 1 hour. 
%       At the end, you should see three figures, each containing a sphere. 
%       The colors on each subsequent sphere is sharper. The colors correspond to sulcal depth.
%       Red => sulci and Blue => gyri.
%       Intermediate subject warps are saved as ".mat" files 
%       Final subject warps are saved as "*h.sphere.SD.reg" in the "surf" directory of each subject. 
%   
% This example function 
%   1)  Go through all the subjects specified in the subject_cell, assumes they are in the SUBJECTS_DIR, 
%   2)  Read in surfaces (assumed to be *h.sphere in the default settings), 
%   3)  Compute an atlas, consisting of the mean and variance of cortical geometry, defined by 
%               (a) curvature of inflated surface (*h.inflated.H)
%               (b) sulcal depth (average convexity, *h.sulc)
%               (c) curvature of original surface (*h.curv)
%   4)  If DISPLAY_ATLAS == 1, and by default it is, function will draw the atlas, 
%       which will look like a sphere with mean sulcal depth (average convexity) painted on it.
%   5)  Register the surfaces to the atlas
%   6)  Save warp (.mat file) into WORK_DIR (by default "SD") of each subject directory. 
%   7)  Recompute the atlas based on previous warps
%   8)  If DISPLAY_ATLAS == 1, and by default it is, function will draw the atlas, 
%       which will look like a sphere with mean sulcal depth (average convexity) painted on it.
%       The atlases get increasingly sharper (more defined).
%   9)  Repeat the registration (without initialization from previous warps)
%   10) Save warp (.mat file) into WORK_DIR (by default "SD") of each subject directory.
%   11) Registration is multi-resolution:
%               (a) First  level: uses curvature of inflated surface (*h.inflated.H)
%               (b) Second level: uses sulcal depth (average convexity, *h.sulc)
%               (c) Third  level: uses sulcal depth (average convexity, *h.sulc)
%               (d) Fourth level: uses curvature of original surface (*h.curv)
%   12) Default uses 2 rounds of registration. In our experiments, sometimes 3 is better depending on the data set. 
%       Too many rounds result in overfitting, similar to using too little regularization. 
%       See following papers for more discussion of overfitting effects (http://yeoyeo02.googlepages.com/publications)
%
%           (a) Effects of Registration Regularization and Atlas Sharpness on Segmentation Accuracy.
%               B.T.T. Yeo, M. Sabuncu, R. Desikan, B. Fischl, P. Golland 
%               Medical Image Analysis, In Press, 2008
%
%           (b) Effects of Registration Regularization and Atlas Sharpness on Segmentation Accuracy.
%               B.T.T. Yeo, M. Sabuncu, R. Desikan, B. Fischl, P. Golland
%               Proceedings of the International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 
%               volume 4791 of LNCS, 683--691, 2007
%
%   13) Final subject warps are saved as "*h.sphere.SD.reg" in the "surf" directory of each subject. 
%
%  The Spherical Demons algorithm is based on the following papers (http://yeoyeo02.googlepages.com/publications):
%
%           (a) Spherical Demons: Fast Diffeomorphic Landmark-Free Surface Registration.
%               B.T.T. Yeo, M. Sabuncu, T. Vercauteren, N. Ayache, B. Fischl, P. Golland
%               In preparation, 2008
%
%           (b) Spherical Demons: Fast Surface Registration.
%               B.T.T. Yeo, M. Sabuncu, T. Vercauteren, N. Ayache, B. Fischl, P. Golland
%               Proceedings of the International Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI), 
%               volume 5241 of LNCS, 745--753, 2008
%
%
%   *************************** Important Note ****************************
%   As mentioned before, the displated atlases will look sharper and sharper after each round. 
%   With our current settings, if you look at the individual subject, the warp is actually not big: average 4mm. 
%   However, in our experiments (refer to papers), this leads to the best segmentation and Brodmann area alignment results. 
%   See point (12) above on overfitting. 
%   For different data sets, the parameters might need to be changed a little.
%   
%   If you want to see dramatic warps, you can drastically decrease reg_parms.smooth_displacement_iter, 
%   for example, set reg_parms.smooth_displacement_iter = [2 2 2 2]. 
%   This decrease the amount of smoothing of the final displacements allowing for larger displacements.
%   
%   In such a scenario, the warp will become rather large and will lead to folding and invertibility issues.
%   This is an unfortunate numerical issue despite the theories of diffeomorphisms we use. 
%   While there's inbuilt code to unfold the warps. This is not only slow, and we believe the resulting registration is not as good. 
%   In this case, we recommend 
%
%       (a) Set reg_parms.smooth_velocity = 1. Smoothing the velocity update allows for a smoother update and will eliminate folding.
%
%       (b) Possibly increasing reg_parms.smooth_velocity_iter. This will result in even smoother updates.
%           Current settings of reg_parms.smooth_velocity_iter leads to no folding even with reg_parms.smooth_displacement_iter = [2 2 2 2]
%
%       (c) Increase reg_parms.iter appropriately. With smoother updates, the final displacement can still be as large 
%           (since final displacement is controlled by smooth_displacement_iter), but the same amount of deformation requires larger number of iterations. 
%
%   It is totally possible that a combination of decreasing reg_parms.smooth_displacement_iter and 
%   increasing reg_parms.smooth_velocity_iter can lead to better results than our default settings 
%   because we did not strongly explore the tradeoff in our Spherical Demons paper.
%

function CoregisterSurfaces(hemi, SUBJECTS_DIR, subject_cell, DISPLAY_ATLAS)

% =========================================================================
% Create a default Atlas parameter structure. Feel free to change the parameters we set here 
% or pass them in through the input arguments of the function.
% For more advanced settings, you can go into CreateDefaultAtlasParms to change stuff.
% =========================================================================
if(nargin < 4)
    SD_atlas = CreateDefaultAtlasParm('lh');                %left hemisphere ('lh') or right hemisphere ('rh')
    SD_atlas.parms.SUBJECTS_DIR = fullfile('..','..','example_surfaces'); % directory that contains subjects
    SD_atlas.parms.subject_cell = {'OAS1_0001_MR1' ,...     % a cell of subjects in SUBJECTS_DIR
        'OAS1_0003_MR1' ,...
        'OAS1_0004_MR1' ,...
        'OAS1_0005_MR1' ,...
        'OAS1_0006_MR1'};
    DISPLAY_ATLAS = 1;                                      % DISPLAY_ATLAS displays the atlases created during the co-registration.
                                                            % The atlases will get sharper and sharper as co-registration continues.
else
    SD_atlas = CreateDefaultAtlasParm(hemi);        %left hemisphere ('lh') or right hemisphere ('rh')
    SD_atlas.parms.SUBJECTS_DIR = SUBJECTS_DIR;     % directory that contains subjects
    SD_atlas.parms.subject_cell = subject_cell;     % a cell of subjects in SUBJECTS_DIR
end

% =========================================================================
% Create a default registration parameter structure. Feel free to change the parameters we set here.
% For more advanced settings, you can go into CreateDefaultRegParms to change stuff.
% Note that currently, we have 4 multiresolution levels, 
%               to decrease it, set some of reg_parms.iter values to 0, e.g. reg_parms.iter = [0 15 15 15] will skip the coarsest level. 
%               to increase it, it's more complicated and not recommended... 
% =========================================================================
reg_parms = CreateDefaultRegParms;                   
reg_parms.iter = [15 15 15 15];                     % How many iterations for each multi-resolution level (currently assumed to be 4 levels)

reg_parms.smooth_displacement_iter = [10 10 10 10]; % How many iterations to smooth the FINAL displacement field in the Spherical Demons second step. 
                                                    % More iterations, smoother FINAL warps (acts like elastic prior).

reg_parms.smooth_velocity = 0;                      % Default 0. If set to 1 => smooth velocity UPDATE in Spherical Demons step 1.
reg_parms.smooth_velocity_iter = [4 4 4 4];         % How many iterations to smooth the velocity UPDATE in the Spherical Demons first step.
                                                    % More iterations, smoother UPDATE (acts like fluid prior). 
                                                    % Currently not being used because reg_parms.smooth_velocity set to 0.
                                                    % If you decide to turn on smooth velocity, you might need to
                                                    % increase reg_parms.iter to get the same amount of warp as when smooth velocity is turned off. In our
                                                    % experiments, setting reg_parms.iter = [20 20 20 20] (when smooth_velocity_iter = [4 4 4 4]) 
                                                    % results in a lower harmonic energy AND better objective function value, BUT slower
                                                    % registration of course!
                                                    
reg_parms.final_unfold = 0;                         % Default 0. The FreeSurfer meshes we use are EXTREMELY irregular
                                                    % We get mesh edge length of as small as 1e-7, average length is 1
                                                    % so we start with edges that are almost 0 (we use floating point precision to save space) 
                                                    % As a result to perform the registration, for each multi-resolution
                                                    % level, we interpolate the subject data onto a regular mesh. At
                                                    % the end of finest level (registration is completed), we interpolate the final
                                                    % warp back onto the original FreeSurfer mesh. Unfortunately, this results in
                                                    % small numbers of folded triangles, typically 10-20 faces out of more than 200k
                                                    % faces. Because of the irregularity of FreeSurfer mesh unfolding these triangles 
                                                    % is really hard, so we do NOT try. If you have meshes that are relative regular,
                                                    % you can set reg_parms.final_unfold to 1, so that after interpolating the
                                                    % warps onto the original mesh, it unfolds any folded triangles
                                                    % introduced in the final warp interpolation.
                                                                                                        
num_iterations_coregistration = 3;                  % number of rounds of co-registration. In our experiments, 2 or 3 works ok depending on the dataset.
                                                    
                                                  
% =========================================================================
% Real code begins!                                                    
% =========================================================================                                                                                                        

% Read in uniform meshes.
disp('Read Uniform Meshes');
reg_parms.meshes{1} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic4.tri');
reg_parms.meshes{2} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic5.tri');
reg_parms.meshes{3} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic6.tri');
reg_parms.meshes{4} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic7.tri');

for iter = 1:num_iterations_coregistration
    
    % Build atlas. Note that the default reads in lh.sphere or rh.sphere (depending on the hemisphere)
    disp(['Iteration ' num2str(iter)]);
    disp('=======================================');
    disp(['Normalizing warps ' num2str(iter)]); 
    SD_atlas = SD_NormalizeAtlasWarps(SD_atlas);
    disp(['Computing atlas ' num2str(iter)]); 
    SD_atlas = SD_CreateAtlasFromRegisteredSurfaces(SD_atlas);
    
    %Display Intermediate atlas
    if(DISPLAY_ATLAS) 
        TrisurfMeshData(SD_atlas.basic_atlas{3}, SD_atlas.basic_atlas{3}.mean, 1); shading interp; title(['Atlas Iteration ' num2str(iter) ': mean "sulcal depth"']);
        view(250, 20);
    end
    
    for subject = 1:length(SD_atlas.parms.subject_cell)
        
        % Read subject mesh
        disp(['Read subject ' num2str(subject)]);
        SD_atlas.parms.SUBJECT = SD_atlas.parms.subject_cell{subject};
        if(~isfield(SD_atlas.parms, 'read_surface'))
            SD_sbjMesh = MARS2_readSbjMesh(SD_atlas.parms);
        else
            SD_sbjMesh = feval(SD_atlas.parms.read_surface, SD_atlas.parms);
        end
    
        % Register subject
        disp(['Register subject ' num2str(subject) ' to atlas']);
        tic
        reg_parms.sbjWarp = []; %no initialization
        reg_parms = SD_registerAtlas2SphereMultiRes(SD_sbjMesh, SD_atlas, reg_parms);
        toc
        
        % Save warp
        sbjWarp = reg_parms.sbjWarp;
        warp_filename = [SD_atlas.parms.hemi '.' num2str(iter) '.warp.mat'];
        save_path = fullfile(SD_atlas.parms.SUBJECTS_DIR, SD_atlas.parms.SUBJECT , SD_atlas.parms.WORK_DIR, warp_filename);
        
        disp(['Saving warp in ' save_path]);
        [not_used1, not_used2] = system(['mkdir ' fullfile(SD_atlas.parms.SUBJECTS_DIR, SD_atlas.parms.SUBJECT , SD_atlas.parms.WORK_DIR)]);
        save(save_path, 'sbjWarp');
        
        % Save warp if final iteration into surf directory of subject.
        if(iter == num_iterations_coregistration)
            final_warp_filename = [SD_atlas.parms.hemi '.sphere.SD.reg'];
            save_path = fullfile(SD_atlas.parms.SUBJECTS_DIR, SD_atlas.parms.SUBJECT , 'surf', final_warp_filename);
            disp(['Saving final warp in ' save_path]);
            if(SD_atlas.parms.flipFacesBool)
                write_surf(sbjWarp.curr_vertices', int32(transpose([SD_sbjMesh.faces(3, :); SD_sbjMesh.faces(2, :); SD_sbjMesh.faces(1, :)]) - 1), save_path);
            else
                write_surf(sbjWarp.curr_vertices', int32(SD_sbjMesh.faces' - 1), save_path);
            end
        end
    end
    
    % set the warp filename so that in the next iteration of creating the
    % atlas, SD_CreateAtlasFromRegisteredSurfaces will use previous warps.
    SD_atlas.parms.warp_filename = warp_filename;
end

%Display final atlas
if(DISPLAY_ATLAS) 
    disp(['Normalizing warps ' num2str(iter)]); 
    SD_atlas = SD_NormalizeAtlasWarps(SD_atlas);
    disp('Computing Final Atlas for display purposes (not used for registration.)');
    SD_atlas = SD_CreateAtlasFromRegisteredSurfaces(SD_atlas);
    TrisurfMeshData(SD_atlas.basic_atlas{3}, SD_atlas.basic_atlas{3}.mean, 1); shading interp; title(['Atlas Final Iteration: mean "sulcal depth"']);
    view(250, 20);
end






