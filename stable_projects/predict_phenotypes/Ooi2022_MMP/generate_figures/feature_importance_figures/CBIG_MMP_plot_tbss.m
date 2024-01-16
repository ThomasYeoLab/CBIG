function CBIG_MMP_plot_tbss(importance, save_dir, prefix, skeleton_mask)

% CBIG_MMP_plot_tbss(importance, save_dir, prefix, skeleton_mask)
% Creates feature importance plots on tbss skeleton. 
%
% Input:
% - importance
% A #tbss_voxels-length vector with importance values.
%
% - save_dir
% The directory in which the images should be saved.
%
% - prefix
% The prefix to append to the saved images.
%
% - skeleton_mask
% Binary tbss skeleton mask derived from thresholding FA.
%
% Output: 
% - Images with file name "prefix_tbss.png"
% One image will be generated representing feature importance values.
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Make sure that there are no NaN, skip if detected
if any(isnan(importance))
    fprintf('NaN detected, unable to plot feature importance. \n')
    return
end

% initialize atlases and tbss mask
std_brain = fullfile(getenv('FSLDIR') , 'data', 'standard', 'MNI152_T1_1mm_brain.nii.gz'); 
cbar_file = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
   'Ooi2022_MMP', 'generate_figures', 'feature_importance_figures', 'colorbar_horz.png');
jhu_atlas = fullfile(getenv('FSLDIR'), 'data', 'atlases', 'JHU', ...
    'JHU-ICBM-labels-1mm.nii.gz');
jhu = MRIread(jhu_atlas);
tbss_mask = MRIread(skeleton_mask);
mask_vol = tbss_mask.vol;
% replace mask with importance values
importance_std_norm = importance ./ std(importance); % normalize by std
new_vol = zeros(size(mask_vol));
importance_vol(mask_vol == 1) = importance_std_norm;
% replace atlas tracts with feature importance values
all_loadings = zeros(48,2);
for n = 1:48
    mean_loading = mean(importance_vol(jhu.vol == n & mask_vol == 1));
    new_vol(jhu.vol == n) = mean_loading;
    all_loadings(n,1)  = n;
    all_loadings(n,2)  = mean_loading;
end
% save unsorted loadings
save(fullfile(save_dir, strcat(prefix,'_fi_unsorted.mat')),'all_loadings');
% save loadings
[b,i] = sort(abs(all_loadings(:,2)), 'descend');
csvwrite(fullfile(save_dir, strcat(prefix,'_tbss_feature_imp.csv')),all_loadings(i,:)); 
tbss_mask.vol = new_vol;
tbss_mask.vol(tbss_mask.vol < 0 ) = 0;
MRIwrite(tbss_mask, fullfile(save_dir, strcat(prefix,'_tbss_JHU_positive.nii.gz')));
tbss_mask.vol = new_vol;
tbss_mask.vol(tbss_mask.vol > 0 ) = 0;
tbss_mask.vol = -tbss_mask.vol; % flip values for visualization
MRIwrite(tbss_mask, fullfile(save_dir, strcat(prefix,'_tbss_JHU_negative.nii.gz')));

% generate figure for JHU atlas
output_file = fullfile(save_dir, strcat(prefix,'_tbss_JHU.png'));
cmd = [ 'fsleyes render -of ' output_file ...
    ' -slightbox -zx Z -ss 4 -zr 35 175 -hc -sz 2000 2000 ', std_brain, ' ' ... 
    fullfile(save_dir, strcat(prefix,'_tbss_JHU_positive.nii.gz')) ...
    ' -cm Red-Yellow -dr 0 2 ' ... 
    fullfile(save_dir, strcat(prefix,'_tbss_JHU_negative.nii.gz')) ...
    ' -cm Blue-LightBlue -dr 0 2']; 
system(cmd);
system(CBIG_ReplaceBlackByWhiteBkgCommand(output_file, output_file));
cmd = ['convert -background white ', output_file, ...
    ' ', cbar_file, ' -gravity Center -append ', output_file]; 
system(cmd);

end