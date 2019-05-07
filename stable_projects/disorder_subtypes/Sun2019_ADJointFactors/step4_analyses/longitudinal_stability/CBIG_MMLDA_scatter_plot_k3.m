function CBIG_MMLDA_scatter_plot_k3(prob_x, prob_y, colors, title_name, xlabel_name, ylabel_name, out_name)
% CBIG_MMLDA_scatter_plot_k3(prob_x, prob_y, colors, title_name, xlabel_name, ylabel_name, out_name)
%
% Scatter plot (1x3 subplots) for factor stability. This will generate Figure 6 in the paper.
%
% Input:
%   - prob_x        : N x 3 matrix. N is number of subjects. Factor loadings for x axis.
%   - prob_y        : N x 3 matrix. N is number of subjects. Factor loadings for y axis.
%   - colors        : N x 1 char. Color for each dot. Each row is 'r' or 'g' or others
%   - title_name    : 1 x 3 cell array. Title for each subplot.
%   - xlabel_name   : string. xlabel name.
%   - ylabel_name   : string. ylabel name.
%   - out_name      : output name of the figure without the suffix '.png'
%
% Example:
%   CBIG_MMLDA_scatter_plot_k3(prob_x, prob_y, colors, {'Memory', 'Language', 'EF'}, xlabel_name, ylabel_name, out_name)
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% set axes and text font
set(0, 'DefaultAxesTickDir', 'out');
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontname', 'Arial'); % Times New Roman
set(0, 'DefaultTextFontSize', 16);

f = figure('Position', [100, 600, 1000, 300]);
for idx = 1:3
    subplot(1, 3, idx);
    oldEst = prob_x(:, idx);
    newEst = prob_y(:, idx);
    % Linear fit
    x = oldEst;
    y = newEst;
    C = corrcoef(x, y);
    text(0.6, 0.1, sprintf('$r=%.2f$',C(1, 2)), 'Interpreter', 'latex', 'FontSize', 20);
    text(0.65, 0.95, '$y=x$', 'Interpreter', 'latex', 'FontSize', 20);
    % t-test between old and new
    [~, p] = ttest(oldEst, newEst);
    fprintf('%s -- Old mean: %f; new mean: %f; p = %e\n', title_name{idx}, mean(oldEst), mean(newEst), p);
    % Plot
    hold on;
    % Need to plot separately for each color. If plotting all together with
    % scatter(x, y, dotSize, colors, 'filled'), saved figure is not vector
    colors_unique = unique(colors);
    for idx2 = 1:numel(colors_unique)
        c = colors_unique(idx2);
        ind = colors==c;
        scatter(x(ind), y(ind), 50, c, 'filled');
    end
    plot(0:1, 0:1, 'k', 'LineWidth', 2); % y = x line
    hold off;
    title(title_name{idx}, 'fontweight', 'bold', 'fontsize', 20);
    xlabel(xlabel_name, 'FontSize', 20);
    ylabel(ylabel_name, 'FontSize', 20);
    xlim([0, 1]);
    ylim([0, 1]);
    axis square;
    ax = gca;
    set(ax, 'XTick', [0 0.5 1], 'FontSize', 20);
    set(ax, 'YTick', [0 0.5 1], 'FontSize', 20);
    set(ax, 'TickDir', 'out');
end

hgexport(f, [out_name '.eps']);
eps2xxx([out_name '.eps'], {'png'});