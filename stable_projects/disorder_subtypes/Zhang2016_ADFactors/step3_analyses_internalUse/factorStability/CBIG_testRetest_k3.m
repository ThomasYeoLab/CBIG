function CBIG_testRetest_k3(rid_prob, rid_prob_new, rid_colors, dotSize, lineColor)

% CBIG_testRetest_k3(rid_prob, rid_prob_new, rid_colors, dotSize, lineColor)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

inputname(2)


% Match
prob_old = [];
prob_new = [];
colors = [];
for idx = 1:size(rid_prob_new, 1)
    rid = rid_prob_new(idx, 1);
    % 
    prob_new = [prob_new; rid_prob_new(idx, 2:end)];
    % 
    ridIdx = rid_prob(:, 1)==rid;
    assert(sum(ridIdx)==1); % a match
    prob_old = [prob_old; rid_prob(ridIdx, 2:end)];
    %
    ridIdx = cell2mat(rid_colors(:, 1))==rid;
    assert(sum(ridIdx)==1); % a match
    colors = [colors; rid_colors{ridIdx, 2}];
end
fprintf('N = %d\n', size(prob_old, 1));

% What is the subtype order?
% Offset by 1, the RID column
fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');
order = [2 3 1];
subtypeNames = {'Temporal', 'Subcortical', 'Cortical'};

CBIG_plotSetup;
f = figure('Position', [100, 600, 1000, 300]);
for idx = 1:3
    subplot(1, 3, idx);
    oldEst = prob_old(:, order(idx));
    newEst = prob_new(:, order(idx));
    % Linear fit
    x = oldEst;
    y = newEst;
    C = corrcoef(x, y);
    text(0.6, 0.1, sprintf('$r=%.2f$',C(1, 2)), 'Interpreter', 'latex', 'FontSize', 20);
    text(0.65, 0.95, '$y=x$', 'Interpreter', 'latex', 'FontSize', 20);
    % t-test between old and new
    [~, p] = ttest(oldEst, newEst);
    fprintf('%s -- Old mean: %f; new mean: %f; p = %e\n', subtypeNames{idx}, mean(oldEst), mean(newEst), p);
    % Plot
    hold on;
    % Need to plot separately for each color. If plotting all together with
    % scatter(x, y, dotSize, colors, 'filled'), saved figure is not vector
    colors_unique = unique(colors);
    for idx2 = 1:numel(colors_unique)
        c = colors_unique(idx2);
        ind = colors==c;
        scatter(x(ind), y(ind), dotSize, c, 'filled');
    end
    plot(0:1, 0:1, lineColor, 'LineWidth', 2); % y = x line
    hold off;
    title([subtypeNames{idx} ' Factor'], 'fontweight', 'bold', 'fontsize', 20);
    xlabel('Probability at Baseline', 'FontSize', 20);
    ylabel('Probability in Month 24', 'FontSize', 20);
    xlim([0, 1]);
    ylim([0, 1]);
    axis square;
    ax = gca;
    set(ax, 'XTick', [0 0.5 1], 'FontSize', 20);
    set(ax, 'YTick', [0 0.5 1], 'FontSize', 20);
    set(ax, 'TickDir', 'out');
end

hgexport(f, [inputname(2) '.eps']);