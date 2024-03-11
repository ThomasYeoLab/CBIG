function CBIG_pFIC_plot_figure_6()

% CBIG_pFIC_plot_figure_6()
% This function generates figures shown in Figure 6 of the manuscript
% 
% Input:
%   - None
% Output:
%   - Different panels of Figure 6
% 
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% Figure 6 (A), (B)
high_performance_age = dlmread('../../replication/GUSTO/input/high_performance/validation_age_list.txt');
low_performance_age = dlmread('../../replication/GUSTO/input/low_performance/validation_age_list.txt');

high_performance_behavior_score = -dlmread('../../replication/GUSTO/input/high_performance/validation_behav_list.txt');
low_performance_behavior_score = -dlmread('../../replication/GUSTO/input/low_performance/validation_behav_list.txt');

% compute p-value
[~, p_age] = ttest(high_performance_age - low_performance_age);
[~, p_behav] = ttest(high_performance_behavior_score - low_performance_behavior_score);
disp('--------------------------------------')
disp(['[Age] high-performance group vs low-performance group p-value: ' num2str(p_age)])
disp(['[Overall cognition] high-performance group vs low-performance group p-value: ' num2str(p_behav)])
disp('--------------------------------------')

hold on
boxplot([high_performance_age', low_performance_age'], ...
    [0*ones(1, length(high_performance_age)), ones(1, length(low_performance_age))],'Colors', 'k', ...
    'Labels', {'high-performance', 'low-performance'});
set(findobj('LineStyle', '--'), 'LineStyle', '-')
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',2)
ylim([7.1, 8.1])
ylabel('age (years)')
set(gca,'Fontsize',15,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
set(gca, 'ytick', [7.2, 7.6, 8])
set(gca,'position',[0.1,0.1,0.80,0.8])
set(gcf,'Position',[800,100,700,700])
print('../../figure/fig6a_age_boxplot', '-dsvg', '-r0')
hold off
close all

hold on
boxplot([high_performance_behavior_score', low_performance_behavior_score'], ...
    [0*ones(1, length(high_performance_behavior_score)), ones(1, length(low_performance_behavior_score))], ...
    'Colors', 'k', 'Labels', {'high-performance', 'low-performance'});
set(findobj('LineStyle', '--'), 'LineStyle', '-')
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',2)
ylim([-3.5, 3.5])
ylabel('overall cognition')
set(gca,'Fontsize',15,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
set(gca, 'ytick', [-3, 0, 3])
set(gca,'position',[0.1,0.1,0.80,0.8])
set(gcf,'Position',[800,100,700,700])
print('../../figure/fig6b_behavior_score_boxplot', '-dsvg', '-r0')
hold off
close all

%% Figure 6 (C)
% To save some storage space on our GitHub repo, the frame-by-frame time
% series of excitatory firing rate and synaptic gating variables are
% deleted from the repo. Instead, the final E/I ratio estimates are saved
% as .csv files under test folder. If the user wants to generate the E/I
% ratio from the time series data, use the helper function compute_EI()
% [At the bottom of this function]
EI_high = dlmread('../../replication/GUSTO/reference_output/test/high_performance/EI.csv');
EI_low = dlmread('../../replication/GUSTO/reference_output/test/low_performance/EI.csv');

% std_split is the standard deviation of regional E/I ratio across 100
% splits of permutation test. To generate std_split, compile all 100 splits
% of E/I ratio from the permutation test and arrange them into a 100x68
% matrix and call it EI_null. std_split = std(EI_null)';
std_split = csvread('../../replication/GUSTO/reference_output/GUSTO_split_std.csv');
effect_size = (EI_low - EI_high)./std_split;

roi_list = ones(68, 1);
cohen_d_fs5 = CBIG_pFIC_ROI2fs5(effect_size, roi_list);
CBIG_pFIC_draw_surface_maps(cohen_d_fs5, 'GUSTO_EI_cognition_cohen_d', './', 'cognition')
movefile('./GUSTO_EI_cognition_cohen_d_lh_lateral.png', '../../figure/fig6c_GUSTO_EI_cognition_cohen_d_lh_lateral.png')
movefile('./GUSTO_EI_cognition_cohen_d_lh_medial.png', '../../figure/fig6c_GUSTO_EI_cognition_cohen_d_lh_medial.png')
movefile('./GUSTO_EI_cognition_cohen_d_rh_lateral.png', '../../figure/fig6c_GUSTO_EI_cognition_cohen_d_rh_lateral.png')
movefile('./GUSTO_EI_cognition_cohen_d_rh_medial.png', '../../figure/fig6c_GUSTO_EI_cognition_cohen_d_rh_medial.png')

%% Figure 6 (D)
effect_size_72 = [0; effect_size(1:3, 1); 0; effect_size(4:34, 1); 0; effect_size(35:37, 1); 0; effect_size(38:68, 1)];
network_assignment = CBIG_pFIC_ROI2network(effect_size_72);
boxplot(network_assignment(:, 2), network_assignment(:, 1), 'Colors', 'k', 'Positions', 1:0.8:5.8, 'Widths', 0.5)
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',1.5)
set(gca,'Fontsize',10,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
ylim([1.58 4.02])
set(gca,'ytick', [1.6 2.8 4.0])
ylabel('Cohen''s d')
set(gca,'xticklabel',{'Som', 'Vis', 'DA', 'VA', 'Lim', 'Control', 'Default'})
set(gcf,'Position',[800,100,1000,700])
print('../../figure/fig6d_boxplot', '-dsvg', '-r0')
close all

%% Figure 6 (E)
SA_rank = dlmread('SA_rank_desikan.txt');
r_spearman = corr(SA_rank, effect_size, 'type', 'Spearman');
disp('--------------------------------------')
disp(['Spearman correlation between effect size and SA axis rank: ' num2str(r_spearman)])
disp('--------------------------------------')
hold on
ylim([1.4 4.5])
set(gca, 'ytick', [1.6 2.8 4.0])
xlim([0, 68.2])
set(gca, 'xtick', [1, 34, 68])
coefficients = polyfit(SA_rank, effect_size, 1);
xFit = linspace(0.01, 70, 1000);
yFit = polyval(coefficients , xFit);
plot(SA_rank, effect_size, '.', 'Color', [105, 105, 105]/255, 'MarkerSize', 15)
set(gca, 'TickDir','out')
set(gca, 'LineWidth', 2)
xlabel('S-A axis rank')
ylabel('Cohen''s d')
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
set(gca, 'FontSize', 15)
hold off
print('../../figure/fig6e', '-dsvg', '-r0')
close all
end

%% helper function
% You found it!
function [EI_high, EI_low] = compute_EI()
rE_min = 2.7;
rE_max = 3.3;

S_E = load('../../replication/GUSTO/reference_output/test/high_performance/simulation/S_E_1.mat');
SE_high = S_E.S_E;
S_I = load('../../replication/GUSTO/reference_output/test/high_performance/simulation/S_I_1.mat');
SI_high = S_I.S_I;
r_E = load('../../replication/GUSTO/reference_output/test/high_performance/simulation/r_E_1.mat');
r_E = r_E.r_E;
anomaly = find(min(r_E) < rE_min);
anomaly = [anomaly find(max(r_E) > rE_max)];
SE_high(:, anomaly) = [];
SI_high(:, anomaly) = [];
EI_high = nanmean(SE_high, 2)./nanmean(SI_high, 2);

S_E = load('../../replication/GUSTO/reference_output/test/low_performance/simulation/S_E_1.mat');
SE_low = S_E.S_E;
S_I = load('../../replication/GUSTO/reference_output/test/low_performance/simulation/S_I_1.mat');
SI_low = S_I.S_I;
r_E = load('../../replication/GUSTO/reference_output/test/low_performance/simulation/r_E_1.mat');
r_E = r_E.r_E;
anomaly = find(min(r_E) < rE_min);
anomaly = [anomaly find(max(r_E) > rE_max)];
SE_low(:, anomaly) = [];
SI_low(:, anomaly) = [];
EI_low = nanmean(SE_low, 2)./nanmean(SI_low, 2);

end
