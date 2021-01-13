function CBIG_SPGrad_example_wrapper(output_dir)

% CBIG_SPGrad_example_wrapper(output_dir)
%
% This function will run examples sequentially
% 
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

 
if(~exist(output_dir))
    mkdir(output_dir);
else
    rmdir(output_dir, 's');
    mkdir(output_dir);
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

%% Generating input example data
system([CBIG_CODE_DIR ...
'/utilities/matlab/speedup_gradients/examples/CBIG_SPGrad_create_example_input_data.sh ' output_dir]);

sub_FC = '100';
sub_verts = '200';
downsample = '3.2';
sub = 1;

out_dir = fullfile(output_dir, ['sub' num2str(sub)]);
mkdir(output_dir);
lh_fMRI_files = fullfile(output_dir, 'data_list', 'fMRI_list', ['lh_sub' num2str(sub) '.txt']);
rh_fMRI_files = fullfile(output_dir, 'data_list', 'fMRI_list', ['rh_sub' num2str(sub) '.txt']);
censor_files = fullfile(output_dir, 'data_list', 'censor_list', ['sub' num2str(sub) '.txt']);

lh_ind_surf = 'NONE';
rh_ind_surf = 'NONE';
medial_mask = 'NONE';
if(~exist(fullfile(output_dir,['rh_gradient_distance_matrix.npy'])))
    if(~exist(fullfile(output_dir,['gradients_edge_density.dtseries.nii'])))
        CBIG_SPGrad_RSFC_gradients(lh_fMRI_files, rh_fMRI_files, censor_files, lh_ind_surf, rh_ind_surf, ...
        'fsaverage6', medial_mask, sub_FC, sub_verts, out_dir);
    end
    CBIG_SPGrad_generate_gradient_matrix('fsaverage6', medial_mask, downsample, out_dir);
end
cmd = ['sh ' CBIG_CODE_DIR '/utilities/matlab/speedup_gradients/CBIG_SPGrad_diffusion_embedding.sh ' out_dir ' 100' ];
system(cmd); 
CBIG_SPGrad_upsample_embed_matrix('fsaverage6', medial_mask, 100, out_dir);

