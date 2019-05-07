function CBIG_MMLDA_generate_JetColormap
% This function will generate a txt file and a png image for Jet colormap.
%
% The JetColormap.txt file is used as input of freeview. For more info, please 
% refer to freeview or CBIG_MMLDA_view_vol_slice_with_underlay.m for example use.
% The JetColorbar.png image is a colorbar for the Jet colormap.
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% create R G B for color table
no_bins = 64;
RGB = colormap(jet(no_bins));
RGB = ceil(RGB*255);

% write color table
fileID = fopen('./JetColormap.txt', 'w');
fprintf(fileID, '#No. LabelName     R    G    B    A\n');
for idx = 1:no_bins
    fprintf(fileID, '%3d      lv%d      %3d  %3d  %3d    0\n', idx, idx, RGB(idx, 1), RGB(idx, 2), RGB(idx, 3));
end
fclose(fileID);

% generate color bar
C = reshape(RGB, [1, size(RGB)]);
C = uint8(C);
image(C);
axis off
set(gcf, 'Color', 'None');
saveas(gcf, 'JetColorbar', 'png');