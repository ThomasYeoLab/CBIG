function CBIG_pFIC_plot_figure_2()

% CBIG_pFIC_plot_figure_2()
% This function generates figures shown in Figure 2 of the manuscript
% 
% Input:
%   - None
% Output:
%   - Different panels of Figure 2
% 
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%% Empirical FC
FC_test = csvread('../../replication/HCP/input/FC_test.csv');
load('FCD_color_table.mat')
imagesc(FC_test);
colormap(c3)
set(gca, 'CLim', [0.1, 1])
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])
set(gca, 'LineWidth', 2)
set(gca,'xtick',[])
set(gca,'ytick',[])
title(gca, 'empirical FC')
%print('../../figure/emprical_FC', '-dsvg', '-r0')
close all

%% Figure 2 (D)
load('../../replication/HCP/reference_output/test/TC.mat')
FC = corr(TC');
TC = TC';
FC_size = 68*(68-1)/2;
FCD_run = zeros(FC_size, 1200 - 82);
for j = 1:(1200 - 82)
    TC_section = TC(j:j+82, :);
    FC_section = CBIG_self_corr(TC_section);
    FC_vec_section = FC_section(triu(true(size(FC_section, 1)), 1));
    FCD_run(:, j) = FC_vec_section;
end
FCD_run = corr(FCD_run);
imagesc(FCD_run);
colormap(c3)
set(gca, 'CLim', [0.1, 1])
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])
set(gca, 'LineWidth', 2)
set(gca,'xtick',[])
set(gca,'ytick',[])
title(gca, 'simulated FCD')
print('../../figure/fig2d', '-dsvg', '-r0')


%% Simulated FC
imagesc(FC);
colormap(c3)
set(gca, 'CLim', [0.1, 1])
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])
set(gca, 'LineWidth', 2)
set(gca,'xtick',[])
set(gca,'ytick',[])
title(gca, 'simulated FC')
%print('../../figure/fig2c', '-dsvg', '-r0')


%% Figure 2 (B)
close all
hold on
plot(FC_test(triu(true(68), 1)), FC(triu(true(68), 1)), '.', 'Color', [105, 105, 105]/255, 'MarkerSize', 8)
xlim([-0.2, 1])
ylim([-0.2, 1])
set(gca, 'LineWidth', 2)
set(gca,'Fontsize', 10,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
set(gca, 'ytick', [0, 0.5, 1])
set(gca, 'xtick', [0, 0.5, 1])
xlabel('empirical FC')
ylabel('simulated FC')
coefficients = polyfit(FC_test(triu(true(68), 1)), FC(triu(true(68), 1)), 1);
xFit = linspace(-0.2, 1, 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 3); % Plot fitted line.
hold off
print('../../figure/fig2b', '-dsvg', '-r0')

%% Figure 2 (C)
close all
load('../../replication/HCP/input/TC.mat')
TC([1 5 37 41], :) = [];
FC = corr(TC');
TC = TC';
FC_size = 68*(68-1)/2;
FCD_run = zeros(FC_size, 1200 - 82);
for j = 1:(1200 - 82)
    TC_section = TC(j:j+82, :);
    FC_section = CBIG_self_corr(TC_section);
    FC_vec_section = FC_section(triu(true(size(FC_section, 1)), 1));
    FCD_run(:, j) = FC_vec_section;
end
FCD_run = corr(FCD_run);
imagesc(FCD_run)
colormap(c3)
set(gca, 'CLim', [0.1, 1])
set(gca, 'YTickLabel', [])
set(gca, 'XTickLabel', [])
set(gca, 'LineWidth', 2)
set(gca,'xtick',[])
set(gca,'ytick',[])
title(gca, 'empirical FCD')
print('../../figure/fig2c', '-dsvg', '-r0')

%% Figure 2 (E)
cost = csvread('../../replication/HCP/reference_output/test/test_all.csv');
FC_corr_cost_pFIC = cost(11, :);
FC_L1_cost_pFIC = cost(12, :);
FCD_cost_pFIC = cost(13, :);

cost = csvread('../../replication/HCP/reference_output/test_homogeneous/test_all.csv');
FC_corr_cost_homo = cost(11, :);
FC_L1_cost_homo = cost(12, :);
FCD_cost_homo = cost(13, :);

cost = csvread('../../replication/HCP/reference_output/test_myelin_only/test_all.csv');
FC_corr_cost_myelin = cost(11, :);
FC_L1_cost_myelin = cost(12, :);
FCD_cost_myelin = cost(13, :);

cost = csvread('../../replication/HCP/reference_output/test_rsfc_gradient_only/test_all.csv');
FC_corr_cost_rsfc = cost(11, :);
FC_L1_cost_rsfc = cost(12, :);
FCD_cost_rsfc = cost(13, :);

total_cost_pFIC = FC_corr_cost_pFIC + FC_L1_cost_pFIC + FCD_cost_pFIC;
total_cost_homo = FC_corr_cost_homo + FC_L1_cost_homo + FCD_cost_homo;
total_cost_myelin = FC_corr_cost_myelin + FC_L1_cost_myelin + FCD_cost_myelin;
total_cost_rsfc = FC_corr_cost_rsfc + FC_L1_cost_rsfc + FCD_cost_rsfc;

boxplot([total_cost_pFIC' total_cost_rsfc' total_cost_myelin' total_cost_homo'], ...
    'Colors', 'k', 'Labels', {'pFIC', 'RSFC gradient', 'T1w/T2w', 'homogeneous'});
ylim([0.5 2.1])
set(findobj('LineStyle', '--'), 'LineStyle', '-')
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'Fontsize',10,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
set(gca, 'ytick', [0.5 1 1.5 2])
set(gca,'position',[0.1,0.1,0.80,0.8])
set(gcf,'Position',[800,100,700,600])
ylabel('test cost (lower = better)')
print('../../figure/fig2g', '-dsvg', '-r0')
close all

end
