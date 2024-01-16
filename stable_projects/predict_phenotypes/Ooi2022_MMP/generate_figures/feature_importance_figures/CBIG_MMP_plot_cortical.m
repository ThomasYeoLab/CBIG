function CBIG_MMP_plot_cortical(importance, save_dir, prefix)

% CBIG_MMP_plot_cortical(importance, save_dir, prefix, mask)
% Creates feature importance plots on cortical surface. Only works for the Schaefer400_17networks parcellation.
%
% Input:
% - importance
% A 400-length vector with importance values in order of the Schaefer400_17networks order.
%
% - save_dir
% The directory in which the images should be saved.
%
% - prefix
% The prefix to append to the saved images.
%
% Output: 
% - 4 images with file name "<prefix>_<hemisphere>.<view>.png"
% An image will be generated for the lh, rh and medial and lateral views.
%
% - prefix_combined.png
% A final image with the lh and rh, lateral and medial and colorbar
% stitched together.
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Make sure that there are no NaN, skip if detected
if any(isnan(importance))
    fprintf('NaN detected, unable to plot feature importance. \n')
    return
end

% convert feature importance to annot file for plotting
importance_std_norm = importance ./ std(importance); % normalize by std
% clip colors for low values: too dark in freesurfer
importance_std_norm(importance_std_norm < 0.1 & importance_std_norm > 0) = 0.1;
importance_std_norm(importance_std_norm > -0.1 & importance_std_norm < 0) = -0.1;
% create annotation file
CBIG_MMP_create_annot_r_overlay_Schaefer400(importance_std_norm, ...
    save_dir, [-2 2], 'hot_cold')
% directory with surf file for visualization
surf_dir = fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', ...
    'surface', 'fs_LR_164k', 'surf');

% plot figures
% save lh lateral
lhl_output_file = fullfile(save_dir, strcat(prefix, '_lh.lateral.png'));
cmd = ['freeview -f ', surf_dir, '/lh.very_inflated:annot=', ...
    save_dir, '/lh_data.annot:edgethickness=0 \ ' ...
    '-zoom 1.78 -ss ', lhl_output_file ];
system(cmd);
system(CBIG_ReplaceBlackByWhiteBkgCommand(lhl_output_file, lhl_output_file));
% save rh lateral
rhl_output_file = fullfile(save_dir, strcat(prefix, '_rh.lateral.png'));
cmd = ['freeview -f ', surf_dir, '/rh.very_inflated:annot=', ...
    save_dir, '/rh_data.annot:edgethickness=0 \ ' ...
    '-cam Azimuth 180 -zoom 1.78 -ss ', rhl_output_file ];
system(cmd);
system(CBIG_ReplaceBlackByWhiteBkgCommand(rhl_output_file, rhl_output_file));
% save lh medial
lhm_output_file = fullfile(save_dir, strcat(prefix, '_lh.medial.png'));
cmd = ['freeview -f ', surf_dir, '/lh.very_inflated:annot=', ...
    save_dir, '/lh_data.annot:edgethickness=0 \ ' ...
    '-cam Azimuth 180 -zoom 1.78 -ss ', save_dir, '/', prefix,'_lh.medial.png'];
system(cmd);
system(CBIG_ReplaceBlackByWhiteBkgCommand(lhm_output_file, lhm_output_file));
% save rh medial
rhm_output_file = fullfile(save_dir, strcat(prefix, '_rh.medial.png'));
cmd = ['freeview -f ', surf_dir, '/rh.very_inflated:annot=', ...
    save_dir, '/rh_data.annot:edgethickness=0 \ ' ...
    '-zoom 1.78 -ss ', rhm_output_file ];
system(cmd);
system(CBIG_ReplaceBlackByWhiteBkgCommand(rhm_output_file, rhm_output_file));

% merge figures and colorbar
final_output_file = strcat(prefix, '_combined.png');
fig_code_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
   'Ooi2022_MMP', 'generate_figures', 'feature_importance_figures');
cbar_file = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
   'Ooi2022_MMP', 'generate_figures', 'feature_importance_figures', 'colorbar_horz.png');

cmd = ['sh ', fig_code_dir, '/CBIG_MMP_merge_cortical.sh ', lhl_output_file, ' ', ...
    rhl_output_file, ' ', lhm_output_file, ' ', rhm_output_file, ' ', save_dir, ...
    ' ', final_output_file ];
system(cmd);
cmd = ['convert -background white ', fullfile(save_dir, final_output_file), ...
    ' ', cbar_file, ' -gravity Center -append ', fullfile(save_dir, final_output_file)]; 
system(cmd);

% remove generated annot files
delete(fullfile(save_dir,'rh_data.annot'))
delete(fullfile(save_dir,'lh_data.annot'))

end