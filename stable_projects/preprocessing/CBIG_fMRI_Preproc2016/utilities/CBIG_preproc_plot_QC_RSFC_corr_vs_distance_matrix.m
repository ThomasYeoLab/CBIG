function CBIG_preproc_plot_QC_RSFC_corr_vs_distance_matrix( QC_vector, FC_matrix, distance_matrix,...
    y_label, output_name )

% CBIG_preproc_plot_QC_RSFC_corr_vs_distance_matrix( QC_vector, FC_matrix, distance_matrix, y_label, output_name )
%
% Plot QC-RSFC correlation (the correlation between quality control measure
% and resting-state functional connecitivity) versus ROIs to ROIs
% volumetric distance.
% 
% Inputs:
%     - QC_vector:
%       An N x 1 vector contains QC measures of all subjects (e.g. mean
%       framewise displacement), where N is number of subjects.
% 
%     - FC_matrix:
%       An M x M x N matrix of functional connectivity for all subjects,
%       where M is number of ROIs.
% 
%     - distance_matrix:
%       An M x M matrix of inter-ROIs volumetric Euclidean distances.
% 
%     - y_label:
%       The ylabel displayed on the plot. For example: 'FD (mean) - RSFC
%       correlation'.
% 
%     - output_name:
%       The full name of output figure, including the path.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% 

%% Check input variables
if size(QC_vector,2) > 1
    error('Input argument ''QC_vector'' should be a column vector');
end

% check the size of matrices
if(length(QC_vector)~=size(FC_matrix, 3))
    error('Number of subjects does not match.')
end
if(size(FC_matrix, 1)~=size(distance_matrix, 1) || size(FC_matrix, 2)~=size(distance_matrix, 2))
    error('Number of ROIs does not match.')
end

[output_dir, ~, ~] = fileparts(output_name);
if(~exist(output_dir, 'dir'))
    mkdir(output_dir);
end

% extract lower triangular indices
tmp_matrix = ones(size(distance_matrix));
tmp_matrix = tril(tmp_matrix, -1);
tril_ind = find(tmp_matrix == 1);
clear tmp_matrix

FC = reshape(FC_matrix, size(FC_matrix, 1)*size(FC_matrix, 2), size(FC_matrix, 3));
FC = FC(tril_ind, :);
distance = reshape(distance_matrix, size(distance_matrix, 1)*size(distance_matrix, 2), size(distance_matrix, 3));
distance = distance(tril_ind, :);

% compute correlation between QC measure and RSFC
QC_vector = QC_vector - mean(QC_vector);
QC_vector = QC_vector ./ sqrt(sum(QC_vector.^2));
FC = bsxfun(@minus, FC, mean(FC, 2));
FC = bsxfun(@times, FC, 1./sqrt(sum(FC.^2, 2)));
QC_RSFC_corr = FC * QC_vector;

% plotting
fig = figure;
hold on
scatter(distance, QC_RSFC_corr, 'r.')
% Create the fitting line
p = polyfit(distance, QC_RSFC_corr,1);
v = polyval(p, distance);
plot(distance, v, '-b', 'LineWidth', 2);
% draw horizontal line 0
xt = get(gca, 'XTick');
plot(xt, zeros(length(xt), 1), 'k', 'LineWidth', 2);
hold off
set(gcf, 'Position', [0, 0, 400, 320])
set(gca, 'LineWidth', 2, 'FontSize', 9)
set(gca, 'YTick', [-1:0.1:1])
xlabel('Inter-node Euclidean Distances (mm)', 'FontSize', 9);
ylabel(y_label, 'FontSize', 9);
% compute stats
mean_QC_RSFC_corr = mean(abs(QC_RSFC_corr));
abs_median = median(abs(QC_RSFC_corr));
RMS_QC_RSFC_corr = sqrt(mean(QC_RSFC_corr.^2));
[r, p] = corrcoef(distance, QC_RSFC_corr);
r = r(1, 2);
p = p(1, 2);

% adjust ylim: round the ylim to be accurate at one decimal place
yl = get(gca, 'ylim');
yll_mod = mod(yl(1)*10, 1);
ylu_mod = mod(yl(2)*10, 1);
if(yll_mod ~= 0)
    yl(1) = floor(yl(1)*10) / 10;
end
if(ylu_mod ~= 0)
    yl(2) = ceil(yl(2)*10) / 10;
end
set(gca, 'ylim', yl);

% write into the title
title_str = {sprintf('Y-axis RMS: %7f, absolute median: %7f', RMS_QC_RSFC_corr, abs_median); ...
    sprintf('X-Y axes Pearson''s r: %7f, p: %7e', r, p)};
ti = title(title_str);
set(ti, 'FontSize', 9)
set(gca, 'TickDir', 'out', 'box', 'off');

% saving
set(gcf, 'PaperPositionMode', 'auto')
print(fig, output_name, '-dpng')
close(fig)


end