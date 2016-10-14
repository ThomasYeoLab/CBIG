function mris_SD_make_template(hemi, input_surface, atlas_name, varargin)

% mris_SD_make_template(<hemi>, <input_surface>, <atlas_name>, <subject>, <subject> ...)
%
% This function goes into the directory specified by the environmental
% variable SUBJECTS_DIR and read in the <hemi> <input_surface> of the list of <subject>
% and creates an atlas and save it in <atlas_name>. Note that the full path
% should be specified inside <atlas_name>.
%
% If the environmental variable SUBJECTS_DIR does not exist, function
% assumes current directory is the SUBJECTS_DIR.
% Note that the function also assumes uniform_mesh_dir is ../../

if(nargin < 4)
   disp(' ');
   disp('mris_SD_make_template(<hemi>, <input_surface>, <atlas_name>, <subject>, <subject> ...)');
   disp('This function goes into the directory specified by the environmental');
   disp('variable SUBJECTS_DIR and read in the <hemi> <input_surface> of the list of <subject>');
   disp('and creates an atlas and save it in <atlas_name>. Note that the full path');
   disp('should be specified inside <atlas_name>.');
   disp(' ');
   disp('If the environmental variable SUBJECTS_DIR does not exist, function');
   disp('assumes current directory is the SUBJECTS_DIR.');
   disp('Note that the function also assumes uniform_mesh_dir is ../../');
   return;
end

if(strcmp(atlas_name(end-3:end), '.mat') == 0)
   error('atlas_name must end in .mat'); 
end

SD_atlas = CreateDefaultFreeSurferAtlas(hemi); 
SD_atlas.parms.surf_filename = input_surface;
SD_atlas.parms.subject_cell = varargin;
SD_atlas.parms.uniform_mesh_dir = fullfile('..','..');

SD_atlas.parms.SUBJECTS_DIR = getenv('SUBJECTS_DIR');
if(isempty(SD_atlas.parms.SUBJECTS_DIR))
    SD_atlas.parms.SUBJECTS_DIR = pwd;
end

SD_atlas.parms

SD_atlas = SD_CreateAtlasFromRegisteredSurfaces(SD_atlas);
save(atlas_name, 'SD_atlas');