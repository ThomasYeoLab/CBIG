function CBIG_SPGrad_upsample_embed_matrix(mesh,medial_mask,num_component,output_dir)

% This function upsample the downsampled diffusion embedding matrix to the original resolution
%
% Input:
%     - mesh:
%       resolution of surface mesh, e.g. 'fsaverage6', 'fs_LR_32k'
%
%     - medial_mask: (#num_vertices x 1 binary vector or 'NONE')
%       the medial area mask. Set it as 'NONE' unless the medial wall area is defined differently
%       as fsaverage medial wall area from FreeSurfer or fs_LR_32k medial wall area from HCP. 
%       If the data has a different medial area mask, please pass in a #num_vertices x 1 binary
%       vector. The medial area should be denoted as 1, others should be denoted as 0.
%
%     - num_component: (double)
%       the number of components for downsampled diffusion embedding matrix, e.g. 100.
%
%     - output_dir:
%       output directory to save the results. The upsampled diffusion embedding matrix will rewrite
%       the downsampled matrix: 
%       <output_dir>/lh_emb_<num_component>_distance_matrix.mat
%       <output_dir>/rh_emb_<num_component>_distance_matrix.mat
%
% Example:
%   CBIG_SPGrad_upsample_embed_matrix('fsaverage6', 'NONE', 100, out_dir)
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_component = num2str(num_component);
outputtmpdir = [output_dir '/tmp'];

if(strcmp(mesh,'fs_LR_32k'))
    lh_surf = CBIG_read_fslr_surface('lh', mesh, 'sphere', 'medialwall.annot');
    rh_surf = CBIG_read_fslr_surface('rh', mesh, 'sphere', 'medialwall.annot');
    if(strcmp(medial_mask,'NONE'))
        medial_mask = [lh_surf.MARS_label == 1;rh_surf.MARS_label == 1];
    end
elseif(~isempty(strfind(mesh,'fsaverage')))
    lh_surf = CBIG_ReadNCAvgMesh('lh', mesh, 'sphere', 'cortex');
    rh_surf = CBIG_ReadNCAvgMesh('rh', mesh, 'sphere', 'cortex');
    if(strcmp(medial_mask,'NONE'))
        medial_mask = [lh_surf.MARS_label==1 rh_surf.MARS_label==1]';
    end
end
lh_downsample_file = fullfile(output_dir,['lh_emb_' num_component '_distance_matrix.mat']);
rh_downsample_file = fullfile(output_dir,['rh_emb_' num_component '_distance_matrix.mat']);

if(~exist(lh_downsample_file))
    error('Cannot find diffusion embedding matrix for left hemisphere');
else
    lh_emb_down = load(lh_downsample_file);
end
if(~exist(rh_downsample_file))
    error('Cannot find diffusion embedding matrix for right hemisphere');
else
    rh_emb_down = load(rh_downsample_file);
end


lh_down_surf = CBIG_SPGrad_read_surf_mesh([outputtmpdir, '/downsampled.lh.sphere']);
rh_down_surf = CBIG_SPGrad_read_surf_mesh([outputtmpdir, '/downsampled.rh.sphere']);

disp('Start upsampling!')
if(size(lh_emb_down.emb,1) == size(lh_surf.vertices,2))
    disp('Left hemisphere has already been upsampled ... skip ...')
else
    [emb, ~, ~] = MARS_linearInterpolate(lh_surf.vertices, lh_down_surf, lh_emb_down.emb');
    emb(:,medial_mask(1:size(lh_surf.vertices,2))) = 0;
    emb = emb';
    save(lh_downsample_file, 'emb', '-v7.3');
    clear emb
end
if(size(rh_emb_down.emb,1) == size(rh_surf.vertices,2))
    disp('Right hemisphere has already been upsampled ... skip ...')
else
    [emb, ~, ~] = MARS_linearInterpolate(rh_surf.vertices, rh_down_surf, rh_emb_down.emb');
    emb(:,medial_mask(size(rh_surf.vertices,2)+1:end)) = 0;
    emb = emb';
    save(rh_downsample_file, 'emb', '-v7.3');
end


disp('########## Done!') 