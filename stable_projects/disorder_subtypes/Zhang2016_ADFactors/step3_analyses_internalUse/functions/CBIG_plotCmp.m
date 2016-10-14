function f = CBIG_plotCmp(A, xlimCmd)

% f = CBIG_plotCmp(A, xlimCmd)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

no_comparisons = size(A, 1);

y_max = (no_comparisons-1)*2+1+1;

% Decide figure height
if no_comparisons == 3 % K = 2
    dy_fig = 295;
elseif no_comparisons == 9 % K = 3
    dy_fig = 825;
elseif no_comparisons == 18 % K = 4
    dy_fig = 1645;
end

f = figure('Position', [100, 100, 1000, dy_fig]);
hold on;
% Comparisons
for idx = 1:no_comparisons
    y_pos = (no_comparisons-idx)*2+1;
    h = CBIG_herrorbar(A(idx, 1), y_pos, A(idx, 2), 'r');
    set(h, 'linewidth', 8);
    scatter(A(idx, 1), y_pos, 400, 'b', 'filled');
end
% 0 lines
plot([0 0], [0, y_max], '--', 'LineWidth', 8, 'Color', [1 0.8 0.6]);
hold off;

ax = gca;
set(ax, 'ytick', []);
set(ax, 'ycolor', [1 1 1]);
set(ax, 'LineWidth', 2.3);
set(ax, 'FontName', 'Arial');
set(ax, 'FontSize', 40);
ylim([0, y_max]);
if strcmp(xlimCmd, 'auto')
    % Do nothing
elseif strcmp(xlimCmd, 'auto-center')
    % Compute xlim
    xAbsMax = max(abs([A(:, 1)+A(:, 2); A(:, 1)-A(:, 2)]));
    xlim([-1.05*xAbsMax, 1.05*xAbsMax]);
elseif strcmp(xlimCmd, 'manual-center')
    xAbsMax = input('What is the limit of the x-axis?\n');
    xlim([-xAbsMax, xAbsMax]);
else
    error('Unconfigured');
end

hgexport(f, [inputname(1) '.eps']);