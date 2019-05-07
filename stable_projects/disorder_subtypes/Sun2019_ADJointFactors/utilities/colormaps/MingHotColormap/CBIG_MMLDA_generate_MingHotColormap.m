function CBIG_MMLDA_generate_MingHotColormap
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% create R G B for color table
no_bins = 255*3;

R = zeros(no_bins, 1);
R(1:255) = 1:255;
R(256:end) = 255;

G = zeros(no_bins, 1);
G(256:510) = 1:255;
G(511:end) = 255;

B = zeros(no_bins, 1);
B(511:end) = 1:255;

% write color table
fileID = fopen('./MingHotColormap.txt', 'w');
fprintf(fileID, '#No. LabelName     R    G    B    A\n');
for idx = 1:no_bins
    fprintf(fileID, '%3d      lv%d      %3d  %3d  %3d    0\n', idx, idx, R(idx), G(idx), B(idx));
end
fclose(fileID);

% generate color bar
C = zeros(1, numel(R), 3);
C(:, :, 1) = R;
C(:, :, 2) = G;
C(:, :, 3) = B;
C = uint8(C);
image(C);
axis off
set(gcf, 'Color', 'None');
saveas(gcf, 'MingHotColorbar', 'png');