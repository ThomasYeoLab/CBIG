function CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name)
% CBIG_MMLDA_triangle_plot(bary_coor, colors, vertex_label, out_name)
%
% Generate triangle plot and save it.
%
% Input:
%   - bary_coor     : N x 3 matrix, N is number of subjects. Barycentric coordinates.
%   - colors        : N x 1 char array. Each row is a char 'r' or 'g' or 'b' or others.
%   - vertex_label  : 1 x 3 cell array. Each cell contains label of each vertex.
%   - out_name      : output name of plot
%
% Example:
%   CBIG_MMLDA_triangle_plot(bary_coor, colors, {'F1', 'F2', 'F3'}, out_name)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% create triangulation
P = [0 1/sqrt(3); 0.5 -0.5/sqrt(3); -0.5 -0.5/sqrt(3)];
T = [1 2 3]; % connectivity
TR = triangulation(T, P); % create the triangulation

% draw triangle
figure;
hold on;
plot([P(1, 1) P(2, 1)], [P(1, 2) P(2, 2)], 'Color', 'k', 'LineWidth', 2);
plot([P(2, 1) P(3, 1)], [P(2, 2) P(3, 2)], 'Color', 'k', 'LineWidth', 2);
plot([P(3, 1) P(1, 1)], [P(3, 2) P(1, 2)], 'Color', 'k', 'LineWidth', 2);

% vertex labels
text('String', vertex_label{1}, 'FontSize', 20, 'Position', P(1, :));
text('String', vertex_label{2}, 'FontSize', 20, 'Position', P(2, :));
text('String', vertex_label{3}, 'FontSize', 20, 'Position', P(3, :));

% scatter
ti = ones(size(bary_coor, 1), 1);
cart_coor = barycentricToCartesian(TR, ti, bary_coor);

% Need to plot separately for each color. If plotting all together with
% scatter(x, y, dotSize, colors, 'filled'), saved figure is not vector
colors_unique = unique(colors);
for j = 1:numel(colors_unique)
    c = colors_unique(j);
    ind = colors==c;
    scatter(cart_coor(ind, 1), cart_coor(ind, 2), 100, c, 'filled');
end

hold off;
axis off;

hgexport(gcf, out_name);
eps2xxx([out_name '.eps'], {'png'})
