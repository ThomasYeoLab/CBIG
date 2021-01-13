function CBIG_IndCBM_create_surf_gifti(mesh, output_dir, varargin)

% CBIG_IndCBM_create_surf_gifti(mesh, output_dir)
%
% This function will create .func.gii files of different meshes for cifti
% templates creating.
%
% Input:
%
%     - mesh:
%           Mesh name. Can be 'fsaverage5', 'fsaverage6', 'fsaverage', 
%           'fs_LR_32k', 'fs_LR_164k' or other customized mesh. Need to
%           pass in optional input 'lh' and 'rh' when using customized mesh
%           name. 
%
%     - output_dir: 
%           Path of output folder. Output will be saved under this folder
%           as L.<mesh>.func.gii and R.<mesh>.func.gii.
%
% Optional input:
%
%     - lh: (double or string)
%           Vertex number of the left hemisphere. Only need to pass in when 
%           mesh name is not 'fsaverage5', 'fsaverage6', 'fsaverage', 
%           'fs_LR_32k' or 'fs_LR_164k'.
%
%     - rh: (double or string)
%           Vertex number of the right hemisphere. Only need to pass in  
%           when mesh name is not 'fsaverage5', 'fsaverage6', 'fsaverage', 
%           'fs_LR_32k' or 'fs_LR_164k'.
%
% Example:
% CBIG_IndCBM_create_surf_gifti('fsaverage5', output_dir)
% CBIG_IndCBM_create_surf_gifti('individual_surface', output_dir, 'lh', '123456', 'rh', '123457')
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

pnames = {'lh' 'rh'};
dflts =  {0 0};
[N_lh, N_rh] = internal.stats.parseArgs(pnames, dflts, varargin{:});

if(ischar(N_lh))
    N_lh = str2num(N_lh);
end
if(ischar(N_rh))
    N_rh = str2num(N_rh);
end

if(strcmp(mesh, 'fsaverage5'))
    N_lh = 10242;
    N_rh = 10242;
elseif(strcmp(mesh, 'fsaverage6'))
    N_lh = 40962;
    N_rh = 40962;
elseif(strcmp(mesh, 'fsaverage'))
    N_lh = 163842;
    N_rh = 163842;
elseif(strcmp(mesh, 'fs_LR_32k'))
    N_lh = 32492;
    N_rh = 32492;
elseif(strcmp(mesh, 'fs_LR_164k'))
    N_lh = 163842;
    N_rh = 163842;
else
    if(~N_lh || ~N_rh)
        error('Missing number of vertices.');
    end
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
g = gifti(fullfile(CBIG_CODE_DIR, 'data', 'templates', 'surface', 'fs_LR_32k', 'example_func', 'Parcels_L.func.gii'));
metric = single(ones(N_lh, 1));
g.cdata = metric;
save(g, fullfile(output_dir, ['L.' mesh '.func.gii']));
metric = single(ones(N_rh, 1));
g.cdata = metric;
save(g, fullfile(output_dir, ['R.' mesh '.func.gii']));

end