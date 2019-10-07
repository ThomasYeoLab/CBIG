function CBIG_ASDf_plotCmp(A, xlimCmd, isSex, output_name)
% CBIG_ASDf_plotCmp(A, xlimCmd, isSex, output_name)
% 
% This function plots GLM results (estimated difference between factors and
% standard error). If output_name is specified, the plot will be saved.
%
% Input:
%     - A:
%           Outputs from hypothesis test, which consists of estimated 
%           difference between factors, standard error and p-value.
%     - xlimCmd:
%           Options to configure the x-axis. If 'auto' is entered, do
%           nothing; if 'auto-center' is entered, the plot will be
%           auto-centered ; if 'manual-center' is entered, the user will be
%           asked to input the limit of x-axis, and the plot will be
%           centered according to this limit.
%     - isSex:
%           Binary variable to indicate whether it is to plot results for
%           sex comparison. 1 indicates yes, 0 indicates no. Because for
%           sex comparison the x-axis needs to be flipped.
%     - output_name (optional):
%           File name (including full path) that will be used to save the plot
%
% Example:
%       CBIG_ASDf_plotCmp(baseline, 'auto-center', 1, 'sex_plot')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

no_comparisons = size(A, 1);
y_max = (no_comparisons-1)*2+1+1;

if isSex
    A = -A;
end

%% Decide figure height
if no_comparisons == 1 % K = 2
    dy_fig = 400/3;
elseif no_comparisons == 3 % K = 3
    dy_fig = 1000/3;
elseif no_comparisons == 6 % K = 4
    dy_fig = 1645/3;
end

figure('Position', [100, 100, 1000, dy_fig]);
hold on;

%% Comparisons
for idx = 1:no_comparisons
    y_pos = (no_comparisons-idx)*2+1;
    h = herrorbar(A(idx, 1), y_pos, A(idx, 2), 'r');
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

%% Center the plot
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

%% Save the plot
if nargin > 3 && ~isempty(output_name)
    hgexport(gcf, output_name);
    eps2xxx([output_name '.eps'], {'png'});
end


