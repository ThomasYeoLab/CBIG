function CBIG_pFIC_plot_figure_5()

% CBIG_pFIC_plot_figure_5()
% This function generates figures shown in Figure 5 of the manuscript
% 
% Input:
%   - None
% Output:
%   - Different panels of Figure 5
% 
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Figure 5 (A) - (C)
% E/I
% To save some storage space on our GitHub repo, the frame-by-frame time
% series of excitatory firing rate and synaptic gating variables are
% deleted from the repo. Instead, the final E/I ratio estimates are saved
% as .csv files under test folder. If the user wants to generate the E/I
% ratio from the time series data, use the helper function compute_EI()
% [At the bottom of this function]
high_performance_EI = zeros(68, 14);
low_performance_EI = zeros(68, 14);
for i = 1:14
   high_performance_EI(:, i) = dlmread(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/high_performance/test/' num2str(i) '/EI.csv']); 
   low_performance_EI(:, i) = dlmread(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/low_performance/test/' num2str(i) '/EI.csv']); 
end
high_behav_EI_mean = mean(high_performance_EI)';
low_behav_EI_mean = mean(low_performance_EI)';
EI = [high_behav_EI_mean low_behav_EI_mean];
hold on
boxplot(EI,'Colors', 'k', 'Labels', {'high-performance', 'low-performance'});
set(findobj('LineStyle', '--'), 'LineStyle', '-')
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'Fontsize',25,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
for i = 1:14
   plot([EI(i, 1), EI(i,2)], 'o-', 'color', [105, 105, 105]/255, 'LineWidth', 1.5, ... 
       'MarkerEdgeColor', 'r', 'MarkerSize', 5)
end
set(gca,'box','off')
set(gca, 'ytick', [2.0, 2.4, 2.8])
ylabel('mean cortical E/I ratio');
set(gca,'position',[0.1,0.1,0.80,0.8])
set(gcf,'Position',[800,100,850,700])
set(gca, 'FontSize', 15)
hold off
print('../../figure/fig5c_EI_ratio_boxplot', '-dsvg', '-r0')
close all

% Age
high_performance_age = zeros(14, 1);
low_performance_age = zeros(14, 1);
for i = 1:14
   age = dlmread(['../../replication/PNC/cognition_effect/input/' ... 
       'high_performance/' num2str(i) '/validation_subject_age.txt']); 
   high_performance_age(i) = mean(age);
end
for i = 1:14
   age = dlmread(['../../replication/PNC/cognition_effect/input/' ... 
       'low_performance/' num2str(i) '/validation_subject_age.txt']); 
   low_performance_age(i) = mean(age);
end
hold on
boxplot([high_performance_age/12, low_performance_age/12],'Colors', 'k', ...
    'Labels', {'high-performance', 'low-performance'});
set(findobj('LineStyle', '--'), 'LineStyle', '-')
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',2)
ylim([8, 24])
set(gca,'Fontsize',25,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
for i = 1:14
   plot([high_performance_age(i)/12, low_performance_age(i)/12], 'o-', ...
       'color', [105, 105, 105]/255, 'LineWidth', 1.5, 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
end
set(gca,'box','off')
set(gca, 'ytick', [10, 16, 22])
ylabel('age (years)')
set(gca,'position',[0.1,0.1,0.80,0.8])
set(gcf,'Position',[800,100,850,700])
set(gca, 'FontSize', 15)
hold off
print('../../figure/fig5a_age_boxplot', '-dsvg', '-r0')
close all

% Behavior
high_performance_score = zeros(14, 1);
low_performance_score = zeros(14, 1);
for i = 1:14
   score = dlmread(['../../replication/PNC/cognition_effect/input' ...
       '/high_performance/' num2str(i) '/validation_subject_behav.txt']); 
   high_performance_score(i) = mean(score);
end
for i = 1:14
   score = dlmread(['../../replication/PNC/cognition_effect/input' ...
       '/low_performance/' num2str(i) '/validation_subject_behav.txt']); 
   low_performance_score(i) = mean(score);
end
hold on
boxplot([high_performance_score, low_performance_score],'Colors', 'k', ...
    'Labels', {'high-performance', 'low-performance'});
set(findobj('LineStyle', '--'), 'LineStyle', '-')
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',2)
ylim([-2, 2])
set(gca,'Fontsize',25,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
for i = 1:14
   plot([high_performance_score(i), low_performance_score(i)], 'o-', ...
       'color', [105, 105, 105]/255, 'LineWidth', 1.5, 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
end
set(gca,'box','off')
set(gca, 'ytick', [-1.5, 0, 1.5])
ylabel('overall accuracy')
set(gca,'position',[0.1,0.1,0.80,0.8])
set(gcf,'Position',[800,100,850,700])
set(gca, 'FontSize', 15)
hold off
print('../../figure/fig5b_overall_acc_boxplot', '-dsvg', '-r0')
close all

% compute p-value
[~, p_EI] = ttest(high_behav_EI_mean - low_behav_EI_mean);
[~, p_age] = ttest(high_performance_age/12 - low_performance_age/12);
[~, p_score] = ttest(high_performance_score - low_performance_score);
disp('--------------------------------------')
disp(['[E/I ratio] high-performance group vs low-performance group p-value: ' num2str(p_EI)])
disp(['[Age] high-performance group vs low-performance group p-value: ' num2str(p_age)])
disp(['[Overall accuracy] high-performance group vs low-performance group p-value: ' num2str(p_score)])
disp('--------------------------------------')

%% Figure 5 (D)
cohen_d = zeros(68, 1);
t = zeros(68, 1);
for i = 1:68
    mean_cohen_d(i, 1) = mean(low_performance_EI(i, :) - high_performance_EI(i, :));
    std_cohen_d(i, 1) = std(low_performance_EI(i, :) - high_performance_EI(i, :));
end
cohen_d = mean_cohen_d./std_cohen_d;
roi_list = ones(68, 1);
cohen_d_fs5 = CBIG_pFIC_ROI2fs5(cohen_d, roi_list);
CBIG_pFIC_draw_surface_maps(cohen_d_fs5, 'PNC_EI_cognition_cohen_d', './', 'cognition')
movefile('./PNC_EI_cognition_cohen_d_lh_lateral.png', '../../figure/fig5d_PNC_EI_cognition_cohen_d_lh_lateral.png')
movefile('./PNC_EI_cognition_cohen_d_lh_medial.png', '../../figure/fig5d_PNC_EI_cognition_cohen_d_lh_medial.png')
movefile('./PNC_EI_cognition_cohen_d_rh_lateral.png', '../../figure/fig5d_PNC_EI_cognition_cohen_d_rh_lateral.png')
movefile('./PNC_EI_cognition_cohen_d_rh_medial.png', '../../figure/fig5d_PNC_EI_cognition_cohen_d_rh_medial.png')

%% Figure 5 (E)
cohen_d_72 = [0; cohen_d(1:3, 1); 0; cohen_d(4:34, 1); 0; cohen_d(35:37, 1); 0; cohen_d(38:68, 1)];
network_assignment = CBIG_pFIC_ROI2network(cohen_d_72);
boxplot(network_assignment(:, 2), network_assignment(:, 1), 'Colors', 'k', 'Positions', 1:0.8:5.8, 'Widths', 0.5)
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',1.5)
set(gca,'Fontsize',10,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
ylim([0.72 1.38])
set(gca,'ytick', [0.75 1.05 1.35])
ylabel('Cohen''s d')
set(gca,'xticklabel',{'Som', 'Vis', 'DA', 'VA', 'Lim', 'Control', 'Default'})
set(gcf,'Position',[800,100,1000,700])
print('../../figure/fig5e_boxplot', '-dsvg', '-r0')
close all

%% Figure 5 (F)
sa_rank = dlmread('SA_rank_desikan.txt');  
roi_list = ones(68, 1);
sa_rank_fs5 = CBIG_pFIC_ROI2fs5(sa_rank, roi_list);
CBIG_pFIC_draw_surface_maps(sa_rank_fs5/100, 'SA_axis_rank', './', 'cognition')
movefile('./SA_axis_rank_lh_lateral.png', '../../figure/fig5f_SA_axis_rank_lh_lateral.png')
movefile('./SA_axis_rank_lh_medial.png', '../../figure/fig5f_SA_axis_rank_lh_medial.png')
movefile('./SA_axis_rank_rh_lateral.png', '../../figure/fig5f_SA_axis_rank_rh_lateral.png')
movefile('./SA_axis_rank_rh_medial.png', '../../figure/fig5f_SA_axis_rank_rh_medial.png')

%% FIgure 5 (G)
r_spearman = corr(sa_rank, cohen_d, 'type', 'Spearman');
disp('--------------------------------------')
disp(['Spearman correlation between Cohen''s d and SA axis rank: ' num2str(r_spearman)])
disp('--------------------------------------')
hold on
ylim([0.72 1.38])
set(gca, 'ytick', [0.75, 1.05, 1.35])
xlim([0, 68.2])
set(gca, 'xtick', [1, 34, 68])
coefficients = polyfit(sa_rank, cohen_d, 1);
xFit = linspace(0.01, 70, 1000);
yFit = polyval(coefficients , xFit);

plot(sa_rank, cohen_d, '.', 'Color', [105, 105, 105]/255, 'MarkerSize', 15)
set(gca, 'TickDir','out')
set(gca, 'LineWidth', 2)
xlabel('S-A axis rank')
ylabel('Cohen''s d')
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
set(gca, 'FontSize', 15)
hold off
print('../../figure/fig5g', '-dsvg', '-r0')
close all


end

%% helper function
% You found it!
function [high_performance_EI, low_performance_EI] = compute_EI()
rE_min = 2.7;
rE_max = 3.3;
high_performance_EI = zeros(68, 14);
low_performance_EI = zeros(68, 14);

for j = 1:14
    r_E = load(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/high_performance/test/' num2str(j) '/simulation/r_E_1.mat']); r_E = r_E.r_E;
    S_E = load(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/high_performance/test/' num2str(j) '/simulation/S_E_1.mat']); S_E = S_E.S_E;
    S_I = load(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/high_performance/test/' num2str(j) '/simulation/S_I_1.mat']); S_I = S_I.S_I;

    anomaly = find(min(r_E) < rE_min);
    anomaly = [anomaly find(max(r_E) > rE_max)];
    S_E(:, anomaly) = [];
    S_I(:, anomaly) = [];

    SE = nanmean(S_E, 2);
    SI = nanmean(S_I, 2);

    high_performance_EI(:, j) = SE./SI;
end

for j = 1:14
    r_E = load(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/low_performance/test/' num2str(j) '/simulation/r_E_1.mat']); r_E = r_E.r_E;
    S_E = load(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/low_performance/test/' num2str(j) '/simulation/S_E_1.mat']); S_E = S_E.S_E;
    S_I = load(['../../replication/PNC/cognition_effect/reference_output' ... 
        '/low_performance/test/' num2str(j) '/simulation/S_I_1.mat']); S_I = S_I.S_I;

    anomaly = find(min(r_E) < rE_min);
    anomaly = [anomaly find(max(r_E) > rE_max)];
    S_E(:, anomaly) = [];
    S_I(:, anomaly) = [];

    SE = nanmean(S_E, 2);
    SI = nanmean(S_I, 2);
    
    low_performance_EI(:, j) = SE./SI;
end
end