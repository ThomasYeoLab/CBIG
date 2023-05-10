function CBIG_preproc_multiecho_QC_greyplot(before_MEICA, after_MEICA, FDpath, output)

% CBIG_preproc_multiecho_QC_greyplot(before_MEICA, after_MEICA, FDpath, subject_dir, bold)
%
% This function generates greyplot for multi-echo QC. Three sub-plots will
% be generated. 1)greyplot for image before MEICA, 2)greyplot for image
% after MEICA, and 3)FD 
%
% Inputs:
%     - before_MEICA:
%       This is tedana output without ME-ICA (without denoising).
%       The absolute file path of nifti file (full path).
%       This should be found in tedana output folder with exact same name 
%       in the following example.
%       E.g. '<path-to-image>/desc-optcom_bold.nii.gz'
%     - after_MEICA:
%       This is tedana output with ME-ICA (with denoising).
%       The absolute file path of nifti file (full path).
%       This should be found in tedana output folder with exact same name 
%       in the following example.
%       E.g. '<path-to-image>/desc-optcomDenoised_bold.nii.gz'
%     - FDpath:
%       The absolute file path to FDRMS file
%       E.g. '<path-to-image>/sub005_bld001_e1_rest_skip4_stc_motion_outliers_FDRMS'
%     - output:
%       Absolute file path and file name to output, which is the QC greyplot
%       E.g. '<path-to-output>/sub005_bold001_greyplot.png'
%
% Example:
%     CBIG_preproc_multiecho_QC_greyplot(...
%     '<path-to-image>/desc-optcom_bold.nii.gz', ...
%     '<path-to-image>/desc-optcomDenoised_bold.nii.gz', ...
%     '<path-to-FD>/sub005_bld001_e1_rest_skip4_stc_motion_outliers_FDRMS', ...
%     '<path-to-output>/sub005_bold001_greyplot.png')
% 
% Written by Lyu Xingyu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

before_MEICA = MRIread(before_MEICA);
after_MEICA = MRIread(after_MEICA);
before_MEICA = before_MEICA.vol;
after_MEICA = after_MEICA.vol;
before_MEICA = reshape(before_MEICA, ...
    size(before_MEICA,1)*size(before_MEICA,2)*size(before_MEICA,3),size(before_MEICA,4));
after_MEICA = reshape(after_MEICA,size(after_MEICA,1)*size(after_MEICA,2)*size(after_MEICA,3),size(after_MEICA,4));
GS = mean(before_MEICA, 1);
FD = transpose(load(FDpath));

subplot(3,1,1);
xtick_vec = int16(linspace(1, length(FD), 6));
plot(1:length(FD), FD, 'r', 'LineWidth', 2);
ax = subplot(3,1,1);
set(ax, 'Units', 'points', 'YTick', [0 roundn(min((max(FD)/2),0.25),-2) roundn(min(max(FD),0.5),-2)], ...
    'YTickLabel', {'0', roundn(min((max(FD)/2),0.25),-2), roundn(min(max(FD),0.5),-2)}, ...
    'YLim', [0 roundn(min(max(FD),0.5),-2)]);
set(ax, 'FontSize', 18);
set(ax, 'Position', [130 600 600 55]);
set(ax, 'Xtick', xtick_vec, 'XTickLabel', ''); 
xlim([1 length(FD)]);
set(get(gca, 'YLabel'), 'String', '\color{red}FD');
set(get(gca, 'YLabel'), 'FontSize', 18);
set(gca, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'box', 'off');

subplot(3,1,2);
STD_gm = mean(std(before_MEICA, [], 2));
before_MEICA_zscore = bsxfun(@minus, before_MEICA, mean(before_MEICA, 2));
before_MEICA_zscore = bsxfun(@rdivide, before_MEICA_zscore, std(before_MEICA, [], 2));
before_MEICA_zscore = transpose(detrend(before_MEICA_zscore'));
before_MEICA = before_MEICA_zscore * STD_gm;
corr_arr = zeros(1,length(before_MEICA));
for i = 1:size(before_MEICA)
    corr_arr(i) = corr(transpose(GS), transpose(before_MEICA(i, :)));
end
[~, I1] = sort(corr_arr, 'descend');
before_MEICA = before_MEICA(I1, :);
before_MEICA(any(isnan(before_MEICA), 2), :) = [];
imagesc(before_MEICA, [-20 20]); 
col = colorbar;
colormap gray
set(col, 'FontSize', 18);
ax = subplot(3,1,2);
set(ax, 'Units', 'points')
set(ax, 'FontSize', 18);
set(ax, 'Position', [130 325 600 200]);
set(ax, 'Xtick', xtick_vec);
set(ax, 'YTick', [])
xlim([1 length(FD)]);
set(get(ax, 'YLabel'), 'String', 'gm+wm+csf');
set(get(ax, 'YLabel'), 'FontSize', 18);
set(get(ax, 'title'), 'String', ' before MEICA');
set(ax, 'TickDir', 'out');


subplot(3,1,3);
STD_gm = mean(std(after_MEICA, [], 2));
after_MEICA_zscore = bsxfun(@minus, after_MEICA, mean(after_MEICA, 2));
after_MEICA_zscore = bsxfun(@rdivide, after_MEICA_zscore, std(after_MEICA, [], 2));
after_MEICA_zscore = transpose(detrend(after_MEICA_zscore'));
after_MEICA = after_MEICA_zscore * STD_gm;

after_MEICA = after_MEICA(I1, :);
after_MEICA(any(isnan(after_MEICA), 2), :) = [];
imagesc(after_MEICA, [-20 20]); 
col = colorbar;
colormap gray
ay = subplot(3,1,3);
set(ay, 'Units', 'points')
set(ay, 'FontSize', 18);
set(ay, 'Position', [130 50 600 200]);     
set(ay, 'Xtick', xtick_vec);
set(ay, 'YTick', [])
xlim([1 length(FD)]);
set(get(ay, 'YLabel'), 'String', 'gm+wm+csf');
set(get(ay, 'YLabel'), 'FontSize', 18);
set(get(ay, 'title'), 'String', 'after MEICA');
set(ay, 'TickDir', 'out');

set(gcf, 'units', 'inch', 'position', [0 0 12 10]);
set(gcf, 'color', 'w');
Image = getframe(gcf);
imwrite(Image.cdata, fullfile(output));
end

