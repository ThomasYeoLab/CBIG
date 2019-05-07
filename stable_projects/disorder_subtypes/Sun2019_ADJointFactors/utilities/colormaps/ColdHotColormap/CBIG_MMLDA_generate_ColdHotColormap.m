function CBIG_MMLDA_generate_ColdHotColormap
% Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Create R G B for color table 
no_bins = 255*4;

R = zeros(no_bins, 1);
R(1:255) = 0;
R(255+1:255*2) = 0;
R(255*2+1:255*3) = 1:255;
R(255*3+1:end) = 255;

G = zeros(no_bins, 1);
G(1:255) = 255:-1:1;
G(255+1:255*2) = 0;
G(255*2+1:255*3) = 0;
G(255*3+1:end) = 1:255;

B = zeros(no_bins, 1);
B(1:255) = 255;
B(255+1:255*2) = 255:-1:1;
B(255*2+1:255*3) = 0;
B(255*3+1:end) = 0;

% write color table
fileID = fopen('./ColdHotColormap.txt', 'w');
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
set(gcf, 'Color', 'None');
saveas(gcf, 'MingHotColorbar', 'png');
