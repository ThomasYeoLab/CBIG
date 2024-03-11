function CBIG_pFIC_plot_figure_3()

% CBIG_pFIC_plot_figure_3()
% This function generates figures shown in Figure 3 of the manuscript
% 
% Input:
%   - None
% Output:
%   - Different panels of Figure 3
% 
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%% Figure 3 (B)
load('../../replication/Alprazolam/reference_output/test/drug/extrapolation/EI.mat')
EI_drug = EI';
load('../../replication/Alprazolam/reference_output/test/placebo/extrapolation/EI.mat')
EI_placebo = EI';
EI_contrast = EI_placebo - EI_drug;
% read in null distributions obtained from the permutation tests
EI_drug_null = csvread('../../replication/Alprazolam/reference_output/EI_drug_null.csv');
EI_placebo_null = csvread('../../replication/Alprazolam/reference_output/EI_placebo_null.csv');
EI_contrast_null = abs(EI_drug_null - EI_placebo_null); % to compute 2-tail p-value
for i = 1:68
    p(i, 1) = sum(EI_contrast(i,1) < EI_contrast_null(i, :))/100;
    if p(i,1) == 0
        p(i,1) = (0+1)/(100+1);
    end
end
p_fdr = mafdr(p, 'BHFDR', true); % FDR correction
roi_list = p_fdr < 0.05; % 1 means that a ROI survives FDR, 0 means it doesn't
EI_contrast_fs5 = CBIG_pFIC_ROI2fs5(EI_contrast, roi_list);

% draw surface map (excluding ROI(s) failed FDR correction)
CBIG_pFIC_draw_surface_maps(EI_contrast_fs5, 'Alprazolam_EI_contrast', './')
movefile('./Alprazolam_EI_contrast_lh_lateral.png', '../../figure/fig3b_Alprazolam_EI_contrast_lh_lateral.png')
movefile('./Alprazolam_EI_contrast_lh_medial.png', '../../figure/fig3b_Alprazolam_EI_contrast_lh_medial.png')
movefile('./Alprazolam_EI_contrast_rh_lateral.png', '../../figure/fig3b_Alprazolam_EI_contrast_rh_lateral.png')
movefile('./Alprazolam_EI_contrast_rh_medial.png', '../../figure/fig3b_Alprazolam_EI_contrast_rh_medial.png')

% insert the medial wall ROIs
EI_contrast_72 = [0; EI_contrast(1:3); 0; EI_contrast(4:34); 0; EI_contrast(35:37); 0; EI_contrast(38:68)];
network_assignment = CBIG_pFIC_ROI2network(EI_contrast_72);
boxplot(network_assignment(:, 2), network_assignment(:, 1), 'Colors', 'k', 'Positions', 1:0.8:5.8, 'Widths', 0.5)
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',1.5)
set(gca,'Fontsize',10,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
ylim([-0.1 1.3])
set(gca,'ytick', [0, 0.6, 1.2])
ylabel('E/I ratio contrast')
set(gca,'xticklabel',{'Som', 'Vis', 'DA', 'VA', 'Lim', 'Control', 'Default'})
set(gcf,'Position',[800,100,1000,700])
print('../../figure/fig3b_boxplot', '-dsvg', '-r0')
close all

%% Figure 3 (C)
BZR = csvread('../../replication/Alprazolam/input/BZR_density.csv');
BZR_fs5 = CBIG_pFIC_ROI2fs5(BZR, ones(68, 1));
CBIG_pFIC_draw_surface_maps(BZR_fs5, 'BZR_density', './')
movefile('./BZR_density_lh_lateral.png', '../../figure/fig3c_BZR_density_lh_lateral.png')
movefile('./BZR_density_lh_medial.png', '../../figure/fig3c_BZR_density_lh_medial.png')
movefile('./BZR_density_rh_lateral.png', '../../figure/fig3c_BZR_density_rh_lateral.png')
movefile('./BZR_density_rh_medial.png', '../../figure/fig3c_BZR_density_rh_medial.png')

%% Figure 3 (D)
r = corr(EI_contrast, BZR);
disp('--------------------------------------')
disp(['E/I ratio - BZR corr: ' num2str(r)])
disp('--------------------------------------')
hold on
ylim([50, 1050])
set(gca, 'ytick', [100, 550, 1000])
xlim([0.08, 1.25])
set(gca, 'xtick', [0.1, 0.65, 1.2])
coefficients = polyfit(EI_contrast, BZR, 1);
xFit = linspace(0.01, 1.5, 1000);
yFit = polyval(coefficients , xFit);
plot(EI_contrast, BZR, '.', 'Color', [105, 105, 105]/255, 'MarkerSize', 15)
set(gca, 'TickDir','out')
ylabel('BZR density')
xlabel('E/I ratio contrast')
set(gca, 'LineWidth', 2)
plot(xFit, yFit, 'r-', 'LineWidth', 2);
hold off
print('../../figure/fig3d', '-dsvg', '-r0')

close all

end
