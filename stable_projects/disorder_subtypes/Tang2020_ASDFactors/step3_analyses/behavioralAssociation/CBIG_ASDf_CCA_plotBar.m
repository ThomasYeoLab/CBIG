function CBIG_ASDf_CCA_plotBar(score_names, struc_coeff, output_name)
% CBIG_ASDf_CCA_plotBar(score_names, struc_coeff, output_name)
%
% This function plots the CCA behavioral score structure coefficients vs raw behavioral scores as
% bar plot. If output_name is specified, the plot will be saved.
%
% Input:
%     - score_names:
%           Px1 cell array, where P is the number of behavioral scores.
%           Each entry is a string
%     - struc_coeff:
%           Px1 matrix. The CCA structure coefficients of the behavioral scores.
%     - output_name (optional):
%           A string. The file name (including full path) that will be used to save the plot
%
% Example:
%     CBIG_ASDf_CCA_plotBar(score_names, struc_coeff, 'K2F1_RRB_barPlot')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if size(struc_coeff,2) > 1
    error('Input argument ''struc_coeff'' should be a column vector');
end

%% Sort the behavior scores from highest structure coefficient to lowest
[struc_coeff_sorted,I] = sort(struc_coeff,'descend');
score_names_sorted = score_names(I);
P = length(score_names); % number of behavior scores

%% Draw bar plot
figure
set(gcf,'Position',[0, 0, 1500, 700]);
hold on

h = bar(struc_coeff_sorted, 0.5);
set(h, 'FaceColor', 'b', 'EdgeColor', 'b');

set(gca, 'xtick', 1:1:1*P, 'xticklabel', score_names_sorted, 'fontsize', 26)
set(gca, 'TickDir', 'out')

rotateXLabels(gca(), 45);
%set(gca, 'xtick', []);

%% Save the plot
if nargin > 2 && ~isempty(output_name)
    hgexport(gcf, output_name);
    eps2xxx([output_name '.eps'], {'png'});
end

hold off;

end