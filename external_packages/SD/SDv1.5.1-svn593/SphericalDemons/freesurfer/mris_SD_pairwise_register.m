function mris_SD_pairwise_register(fixed_surface, moving_surface, output_surface, disp_iter)

% mris_SD_register(input_surface, moving_surface, output_surface, disp_iter)
% Example call:
%      <When in the same directory as mris_SD_pairwise_register.m>
%      mris_SD_pairwise_register('../../example_surfaces/OAS1_0001_MR1/surf/lh.sphere', '../../example_surfaces/OAS1_0003_MR1/surf/lh.sphere', '../../example_surfaces/OAS1_0001_MR1/surf/lh.sphere.SD.reg');
%   


if(nargin < 3)
   disp(' ');
   disp('mris_SD_pairwise_register(fixed_surface, moving_surface, output_surface, <optional> disp_iter)');
   return;
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Read fixed surface
%%%%%%%%%%%%%%%%%%%%%%%%

% Grab hemi, SUBJECTS_DIR from fixed_surface.
indices = strfind(fixed_surface, filesep);
hemi = fixed_surface(indices(end)+1:indices(end)+2);
SUBJECTS_DIR = fixed_surface(1:indices(end-2)-1);
subject = fixed_surface(indices(end-2)+1:indices(end-1)-1);
surf_dir = fixed_surface(indices(end-1)+1:indices(end)-1); 
if(strcmp(surf_dir, 'surf') == 0)
   error(['mris_SD_register assumes surfaces and data are in the "surf" directory of the subject, not ' surf_dir]); 
end

SD_atlas = CreateDefaultFreeSurferAtlas(hemi);
SD_atlas.parms.uniform_mesh_dir = fullfile('..','..');
SD_atlas.parms.surf_filename = fixed_surface(indices(end)+4:end); 
SD_atlas.parms.SUBJECT = subject;
SD_atlas.parms.SUBJECTS_DIR = SUBJECTS_DIR;

SD_atlas.parms

% Read fixed subject mesh
disp(['Read fixed subject ' num2str(subject) ': ' SD_atlas.parms.hemi '.' SD_atlas.parms.surf_filename]);
disp(['and data: ' SD_atlas.parms.data_filename_cell]);

if(~isfield(SD_atlas.parms, 'read_surface'))
    fixed_sbjMesh = MARS2_readSbjMesh(SD_atlas.parms);
else
    fixed_sbjMesh = feval(SD_atlas.parms.read_surface, SD_atlas.parms);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Read moving surface
%%%%%%%%%%%%%%%%%%%%%%%%

% Grab hemi, SUBJECTS_DIR from moving_surface.
indices = strfind(moving_surface, filesep);
hemi = moving_surface(indices(end)+1:indices(end)+2);
SUBJECTS_DIR = moving_surface(1:indices(end-2)-1);
subject = moving_surface(indices(end-2)+1:indices(end-1)-1);
surf_dir = moving_surface(indices(end-1)+1:indices(end)-1); 
if(strcmp(surf_dir, 'surf') == 0)
   error(['mris_SD_register assumes surfaces and data are in the "surf" directory of the subject, not ' surf_dir]); 
end

% check hemi
if(strcmp(SD_atlas.parms.hemi, hemi) == 0)
   error(['Fixed surface uses hemisphere ' SD_atlas.parms.hemi ', but moving_surface is ' hemi]);
end

SD_atlas.parms.surf_filename = moving_surface(indices(end)+4:end); 
SD_atlas.parms.SUBJECT = subject;
SD_atlas.parms.SUBJECTS_DIR = SUBJECTS_DIR;

SD_atlas.parms

% Read fixed subject mesh
disp(['Read moving subject ' num2str(subject) ': ' SD_atlas.parms.hemi '.' SD_atlas.parms.surf_filename]);
disp(['and data: ' SD_atlas.parms.data_filename_cell]);

if(~isfield(SD_atlas.parms, 'read_surface'))
    moving_sbjMesh = MARS2_readSbjMesh(SD_atlas.parms);
else
    moving_sbjMesh = feval(SD_atlas.parms.read_surface, SD_atlas.parms);
end


% initialize reg_parms (These are the defaults)
reg_parms = CreateDefaultFreeSurferRegParms; 
reg_parms.iter = [15 15 15 15];                     % How many iterations for each multi-resolution level (currently assumed to be 4 levels)

reg_parms.smooth_displacement_iter = [10 10 10 10]; % How many iterations to smooth the FINAL displacement field in the Spherical Demons second step. 
                                                    % More iterations, smoother FINAL warps (acts like elastic prior).

if(nargin == 4)                                                    
    reg_parms.smooth_displacement_iter(:) = str2num(disp_iter);                                                    
end
    
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
                                                    % introduced in the
                                                    % final warp interpolation.


% =========================================================================
% Real code begins!                                                    
% =========================================================================                                                                                                        

% Read in uniform meshes.
disp('Read Uniform Meshes');
reg_parms.meshes{1} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic4.tri');
reg_parms.meshes{2} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic5.tri');
reg_parms.meshes{3} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic6.tri');
reg_parms.meshes{4} = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic7.tri');


% Register subject
disp(['Register moving surface and fixed surface']);
tic
reg_parms.sbjWarp = [] %no initialization
reg_parms = SD_registerPairOfSpheres(fixed_sbjMesh, moving_sbjMesh, SD_atlas, reg_parms);
sbjWarp = reg_parms.sbjWarp;
toc

% save warp
save_path = output_surface;
disp(['Saving final warp in ' save_path]);
if(SD_atlas.parms.flipFacesBool)
    write_surf(sbjWarp.curr_vertices', int32(transpose([fixed_sbjMesh.faces(3, :); fixed_sbjMesh.faces(2, :); fixed_sbjMesh.faces(1, :)]) - 1), save_path);
else
    write_surf(sbjWarp.curr_vertices', int32(fixed_sbjMesh.faces' - 1), save_path);
end



