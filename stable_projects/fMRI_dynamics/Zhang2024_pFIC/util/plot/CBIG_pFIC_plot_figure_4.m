function CBIG_pFIC_plot_figure_4()

% CBIG_pFIC_plot_figure_4()
% This function generates figures shown in Figure 4 of the manuscript
% 
% Input:
%   - None
% Output:
%   - Different panels of Figure 4
% 
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Figure 4 (A)
age = [];
for i = 1:29
    age_training = dlmread(['../../replication/PNC/age_effect/input/' num2str(i) '/training_subject_age.txt']);
    age_validation = dlmread(['../../replication/PNC/age_effect/input/' num2str(i) '/validation_subject_age.txt']);
    age = [age; age_training; age_validation];
end
histogram(age/12, 'LineWidth', 2, 'FaceColor', [105, 105, 105]/255)
set(gca, 'TickDir','out')
set(gca,'box','off')
set(gca, 'LineWidth', 2)
set(gca, 'ytick', [0, 50, 100])
set(gca, 'xtick', [8, 12, 16, 20, 24])
ylabel('number of participants')
xlabel('age (years)')
set(gca, 'FontSize', 15)

print('../../figure/fig4a', '-dsvg', '-r0')
close all

%% Figure 4 (B)
% To save some storage space on our GitHub repo, the frame-by-frame time
% series of excitatory firing rate and synaptic gating variables are
% deleted from the repo. Instead, the final E/I ratio estimates are saved
% as .csv files under test folder. If the user wants to generate the E/I
% ratio from the time series data, use the helper function compute_EI()
% [At the bottom of this function]

age = zeros(29, 1);
EI = zeros(68, 29);
for i = 1:29
    EI(:, i) = csvread(['../../replication/PNC/age_effect/reference_output/test/' num2str(i) ... 
        '/EI.csv']);
    age(i, 1) = mean(dlmread(['../../replication/PNC/age_effect/input/' num2str(i) ...
        '/validation_subject_age.txt']));
end

EI_mean = mean(EI)';
[r_age, p_age] = corr(EI_mean, age);
disp('--------------------------------------')
disp(['E/I vs age corr: ' num2str(r_age)])
disp(['E/I vs age p-value: ' num2str(p_age)])
disp('--------------------------------------')

hold on
ylim([1.8 3])
set(gca, 'ytick', [2 2.4 2.8])
xlim([8 23])
set(gca, 'xtick', [12 16 20])
plot(age/12, EI_mean, '.', 'Color', [105, 105, 105]/255, 'MarkerSize', 15)
mdl = fitlm(age/12, EI_mean, 'linear');
h = plot(mdl);
[~, ci] = predict(mdl, age/12);  
patch([age/12; flipud(age/12)], [ci(:, 1); flipud(ci(:, 2))], ...
    'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
set(h, 'Marker', 'none', 'LineWidth', 2);
legend('off');
title(''); 
set(gca, 'TickDir','out')
set(gca, 'LineWidth', 2)
ylabel('mean cortical E/I ratio')
xlabel('group-average age (years)')
set(gca, 'FontSize', 15)
hold off
print('../../figure/fig4b', '-dsvg', '-r0')
close all

%% Figure 4 (C)
X = age/12;
slope = zeros(68, 2);
for i = 1:68
    lm = fitlm(X, EI(i, :)');
    slope(i, 1) = lm.Coefficients{2, 1};
    slope(i, 2) = lm.Coefficients{2, 4};
end
roi_list = mafdr(slope(:, 2), 'BHFDR', 1) < 0.05;
disp('--------------------------------------')
disp(['Number of significant ROI after FDR correction: ' num2str(sum(roi_list))]);
disp('--------------------------------------')

slope_fs5 = CBIG_pFIC_ROI2fs5(slope(:, 1), roi_list);

% draw surface map (excluding ROI(s) failed FDR correction)
% CBIG_pFIC_draw_surface_maps(slope_fs5, 'PNC_EI_age_slope', './')
% movefile('./PNC_EI_age_slope_lh_lateral.png', '../../figure/fig4c_PNC_EI_age_slope_lh_lateral.png')
% movefile('./PNC_EI_age_slope_lh_medial.png', '../../figure/fig4c_PNC_EI_age_slope_lh_medial.png')
% movefile('./PNC_EI_age_slope_rh_lateral.png', '../../figure/fig4c_PNC_EI_age_slope_rh_lateral.png')
% movefile('./PNC_EI_age_slope_rh_medial.png', '../../figure/fig4c_PNC_EI_age_slope_rh_medial.png')

%% Figure 4 (D)

% insert the medial wall ROIs
slope_72 = [0; slope(1:3, 1); 0; slope(4:34, 1); 0; slope(35:37, 1); 0; slope(38:68, 1)];
network_assignment = CBIG_pFIC_ROI2network(slope_72);
boxplot(network_assignment(:, 2), network_assignment(:, 1), 'Colors', 'k', 'Positions', 1:0.8:5.8, 'Widths', 0.5)
set(findobj(gcf,'tag','Outliers'), 'MarkerEdgeColor', 'black', 'MarkerSize', 5, 'Marker', '+')
set(findobj(gca,'type','line'),'linew',1.5)
set(gca,'Fontsize',10,'TickDir','out','FontWeight','bold')
set(gca,'LineWidth',2)
set(gca,'box','off')
ylim([-0.0455 -0.0265])
set(gca,'ytick', [-0.045 -0.036 -0.027])
ylabel('slope')
set(gca,'xticklabel',{'Som', 'Vis', 'DA', 'VA', 'Lim', 'Control', 'Default'})
set(gcf,'Position',[800,100,1000,700])
print('../../figure/fig4d_boxplot', '-dsvg', '-r0')
close all

%% helper function
% You found it!
function [EI, age] = compute_EI()
    rE_min = 2.7;
    rE_max = 3.3;
    for group_number = 1:29
        r_E = load(['../../replication/PNC/age_effect/reference_output/test/' num2str(group_number) ...
            '/simulation/r_E_1.mat']); r_E = r_E.r_E;
        S_E = load(['../../replication/PNC/age_effect/reference_output/test/' num2str(group_number) ...
            '/simulation/S_E_1.mat']); S_E = S_E.S_E;
        S_I = load(['../../replication/PNC/age_effect/reference_output/test/' num2str(group_number) ...
            '/simulation/S_I_1.mat']); S_I = S_I.S_I;

        age_list = dlmread(['../../replication/PNC/age_effect/input/' num2str(group_number) ...
            '/validation_subject_age.txt']);

        anomaly = find(min(r_E) < rE_min);
        anomaly = [anomaly find(max(r_E) > rE_max)];
        S_E(:, anomaly) = [];
        S_I(:, anomaly) = [];

        SE = nanmean(S_E, 2);
        SI = nanmean(S_I, 2);

        EI(:, group_number) = SE./SI;
        age(group_number, 1) = mean(age_list);
    end
end

end
