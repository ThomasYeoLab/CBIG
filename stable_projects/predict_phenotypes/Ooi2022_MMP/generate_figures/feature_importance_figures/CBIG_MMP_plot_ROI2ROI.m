function CBIG_MMP_plot_ROI2ROI(importance, save_dir, prefix)

% CBIG_MMP_plot_ROI2ROI(importance, save_dir, prefix)
% Creates feature importance plots on network edges. 
%
% Input:
% - importance
% A 79800-length vector with importance values from the 400x400 Schaefer 
% parcellation.
%
% - save_dir
% The directory in which the images should be saved.
%
% - prefix
% The prefix to append to the saved images.
%
% Output: 
% - Images with file name "prefix_roi2roi.png"
% One image will be generated representing feature importance values.
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Make sure that there are no NaN, skip if detected
if any(isnan(importance))
    fprintf('NaN detected, unable to plot feature importance. \n')
    return
end

% prepare feature importance vector for plotting
importance_max_norm = importance ./ std(importance); % normalize by std
imp_mat = CBIG_MMP_FC_vector_2_mat(importance_max_norm);
% plot ROI2ROI matrix 
CBIG_PlotCorrMatNetOrder(400, imp_mat, 'Schaefer_Kong17', [-2 2])
% remove colorbar and text and set size
colorbar('hide')
fig_t = findall(gcf,'Type','text');
delete(fig_t)
set(gca, 'Position', [0 0 1 1]) 
set(gcf,'PaperPosition',[0 0 6 6])
saveas(gcf, fullfile(save_dir, strcat(prefix,'_roi2roi.png'))); close(gcf);

end