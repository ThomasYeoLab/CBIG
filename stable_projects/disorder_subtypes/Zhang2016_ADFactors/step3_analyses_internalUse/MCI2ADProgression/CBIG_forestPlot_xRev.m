function f = CBIG_forestPlot_xRev(A)

% f = CBIG_forestPlot_xRev(A)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

f = figure('Position', [100, 100, 1000, 285]);

hold on;

% 0 lines
plot([0 0], [0, 6], '--', 'LineWidth', 8, 'Color', [1 0.8 0.6]);
% plot(xrange, [0 0], 'Color', [0 0.5 0.5], 'LineWidth', 2);

% Whiskers
h = CBIG_herrorbar(A(3, 1), 1, A(3, 2), 'r');
set(h, 'linewidth', 8);
% Dot
scatter(A(3, 1), 1, 400, 'b', 'filled');

% Whiskers
h = CBIG_herrorbar(A(2, 1), 3, A(2, 2), 'r');
set(h, 'linewidth', 8);
% Dot
scatter(A(2, 1), 3, 400, 'b', 'filled');

% Whiskers
h = CBIG_herrorbar(A(1, 1), 5, A(1, 2), 'r');
set(h, 'linewidth', 8);
% Dot
scatter(A(1, 1), 5, 400, 'b', 'filled');

hold off;

ax = gca;
set(ax, 'ytick', []);
set(ax, 'ycolor', [1 1 1]);
set(ax, 'xtick', -5:1:5);
set(ax, 'LineWidth', 3);
set(ax, 'FontName', 'Arial');
set(ax, 'FontSize', 40);
% x-axis flipped so that left is worse
set(ax, 'xDir', 'Reverse');

xlim([-6 6]);
ylim([0, 6]);

hgexport(f, './conv.eps');