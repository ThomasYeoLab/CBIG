function CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx, output_name)
% CBIG_ASDf_CCA_plotScatter(V, U, R, pVal, factor_idx)
% 
% This function plots scatter plot between factor loading CCA mode and 
% behavioral scores CCA mode. If output_name is specified, the plot will be saved.
%
% Input:
%      - V:
%           Nx1 matrix, where N is the number of subjects. This is the
%           canonical variate (mode) of factor loading.
%      - U:
%           Nx1 matrix. The canonical variate (mode) of behavior scores
%      - R:
%           A scalar. The canonical correlation from CCA
%      - pVal:
%           Permutation test p-value
%      - factor_idx:
%           Index of the current factor
%     - output_name (optional):
%           A string. The file name (including path) that will be used to save the plot
%
% Example:
%     CBIG_ASDf_CCA_plotScatter(V, U, 0.45, 0.004, 1)
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if size(V,2) > 1
    error('Input argument ''V'' should be a column vector');
end

if size(U,2) > 1
    error('Input argument ''U'' should be a column vector');
end

%% Draw scatter plot
figure;
hold on;
scatter(V, U, 'b', 'o', 'filled');
lsline;

y_min = floor(min(U));
y_max = ceil(max(U));
ylim([y_min y_max]);
 
xlabel(['Factor ' num2str(factor_idx) 'CCA Mode'], 'fontsize', 25);
ylabel('Behavior Score CCA Mode', 'fontsize', 25);

%% Add text
txt_r = ['r = ' num2str(R, '%.2f')];
txt_p = ['p = ' num2str(pVal, '%.2e')];

text(0.8, y_min+1, txt_r, 'fontsize', 15);
text(0.8, y_min+0.5, txt_p, 'fontsize', 15);

%% Save the plot
if nargin > 5 && ~isempty(output_name)
    hgexport(gcf, output_name);
    eps2xxx([output_name '.eps'], {'png'});
end

hold off;

