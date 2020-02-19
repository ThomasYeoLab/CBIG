function p_10foldCV = CBIG_MMLDA_10foldCV_wrapper(out_dir)
% p_10foldCV = CBIG_MMLDA_10foldCV_wrapper(out_dir)
%
% Wrapper function to do 10 fold cross validation.
%
% Input:
%   - out_dir   : absolute path to the output directory
%
% Output:
%   - p_10foldCV: correlation p values used for FDR correction later
%
% Example:
%   p_10foldCV = CBIG_MMLDA_10foldCV_wrapper('~/example/')
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

proj_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors'];

% Find the factor's order of each fold
k = 3;
fold_num = 10;
ref_dir = [proj_dir '/step2_MMLDA/results/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1_train1'];
ref_r = CBIG_MMLDA_get_run_num(ref_dir, k);
ref_order = [1 3 2];

avg_corr_all(1) = 1;
order_all(1, :) = ref_order;
for i = 2:fold_num
    in_dir = [proj_dir '/step2_MMLDA/results/visualizeFactors/ADNI2_bl_AD_meanCNstdALL_plus1_train' ...
        num2str(i)];
    in_r = CBIG_MMLDA_get_run_num(in_dir, k);
    [in_order, avg_corr] = CBIG_MMLDA_hungarian_match_factors([ref_dir '/k' num2str(k) '/r' num2str(ref_r)], ...
        [in_dir '/k' num2str(k) '/r' num2str(in_r)], ref_order);
    avg_corr_all(i) = avg_corr;
    order_all(i, :) = in_order;
end

% Stack factor's compositions of all folds together and calculate correlation
prob1_all = [];
prob2_all = [];
cor_all = zeros(fold_num, k);
p_all = zeros(fold_num, k);
for i = 1:fold_num
    % get the atrophy and behavior factor loadings for each test fold 
    rid_file = [proj_dir '/step2_MMLDA/results/BrainBehavior2doc/ADNI2_bl_AD_meanCNstdALL_plus1_RID_test' ...
        num2str(i) '.txt'];
    inf1_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1_train' ...
        num2str(i) '/k3_inf1_ADNI2_bl_AD_meanCNstdALL_plus1_test' num2str(i) '-gamma.dat'];
    rid_prob1_reorder = CBIG_MMLDA_get_factor_loadings(rid_file, inf1_gamma_file, order_all(i, :));
    inf2_gamma_file = [proj_dir '/step2_MMLDA/results/inference/ADNI2_bl_AD_meanCNstdALL_plus1_train' ...
        num2str(i) '/k3_inf2_ADNI2_bl_AD_meanCNstdALL_plus1_test' num2str(i) '-gamma.dat'];
    rid_prob2_reorder = CBIG_MMLDA_get_factor_loadings(rid_file, inf2_gamma_file, order_all(i, :));
    prob1_all = [prob1_all; rid_prob1_reorder(:, 2:end)];
    prob2_all = [prob2_all; rid_prob2_reorder(:, 2:end)];

    % calculate the correlation of each fold and p value
    cor = zeros(1, k);
    pval = zeros(1, k);
    for j = 1:k
        [r, p] = corrcoef(rid_prob1_reorder(:, j+1), rid_prob2_reorder(:, j+1));
        cor(j) = r(1,2);
        pval(j) = p(1,2);
    end
    cor_all(i, :) = cor;
    p_all(i, :) = pval;
end

% calculate the correlation and p value by stacking all folds
[r, p]= corr(prob2_all, prob1_all);
% p_10foldCV = diag(p);
p_10foldCV = p(:);
r_p = [r p];
csvwrite([out_dir '/10foldCV.csv'], r_p);

% 1x3 scatter plot of the correlation for each factor
% set axes and text font
% set(0, 'DefaultAxesTickDir', 'out');
% set(0, 'DefaultAxesFontName', 'Arial');
% set(0, 'DefaultAxesFontSize', 16);
% set(0, 'DefaultTextFontname', 'Arial'); % Times New Roman
% set(0, 'DefaultTextFontSize', 16);
% f = figure('Position', [100, 600, 1000, 300]); % [left bottom width height]
% for factor = 1:k
%     subplot(1, k, factor);
    
%     text(0.6, 0.95, sprintf('$r=%.2f$',r(factor, factor)), 'Interpreter', 'latex', 'FontSize', 20);
%     text(0.5, 0.85, sprintf('$p=%s$',CBIG_MMLDA_pvalue2str(p(factor, factor))), ...
%         'Interpreter', 'latex', 'FontSize', 20);
%     hold on;
%     scatter(prob1_all(:, factor)', prob2_all(:, factor)', 50, 'r', 'filled')
%     para = polyfit(prob1_all(:, factor), prob2_all(:, factor), 1);
%     yfit = polyval(para, prob1_all(:, factor));
%     plot(prob1_all(:, factor), yfit, 'k', 'LineWidth', 2)
%     hold off;
%     xlabel('Pr(Atrophy Factor | Patient)', 'FontSize', 18)
%     ylabel('Pr(Behavioral Factor | Patient)', 'FontSize', 18)
%     if factor == 1
%         title('Episodic Memory', 'fontweight', 'bold', 'fontsize', 20)
%     elseif factor == 2
%         title('Language', 'fontweight', 'bold', 'fontsize', 20)
%     elseif factor == 3
%         title('Executive Function', 'fontweight', 'bold', 'fontsize', 20)
%     end

%     xlim([0 1]);
%     ylim([0 1]);
%     axis square;
%     ax = gca;
%     set(ax, 'XTick', [0 0.5 1], 'FontSize', 20);
%     set(ax, 'YTick', [0 0.5 1], 'FontSize', 20);
%     set(ax, 'TickDir', 'out');
% end
% hgexport(f, [out_dir '/10_fold_CV.eps'])
% eps2xxx([out_dir '/10_fold_CV.eps'], {'png'})

% 3 x 3 scatter plot 
xlabel_set = {'Medial Temporal Loading', 'Lateral Temporal Loading', 'Posterior Cortical Loading'};
ylabel_set = {'Memory Loading', 'Language Loading', 'Executive Loading'};
f = figure('Position', [100, 100, 1000, 1000]); % [left bottom width height]
for row = 1:k
    for col = 1:k
        subplot(k, k, (row-1)*k + col);
        
        text(0.6, 0.95, sprintf('$r=%.2f$',r(row, col)), 'Interpreter', 'latex', 'FontSize', 20);
        text(0.5, 0.85, sprintf('$p=%s$',CBIG_MMLDA_pvalue2str(p(row, col))), 'Interpreter', 'latex', 'FontSize', 20);
        hold on;
        if row == col
            scatter(prob1_all(:, col)', prob2_all(:, row)', 50, 'r', 'filled')
        else
            scatter(prob1_all(:, col)', prob2_all(:, row)', 50, 'b', 'filled')
        end
        para = polyfit(prob1_all(:, col), prob2_all(:, row), 1);
        yfit = polyval(para, prob1_all(:, col));
        plot(prob1_all(:, col), yfit, 'k', 'LineWidth', 2)
        hold off;
        if row == 3
            xlabel(xlabel_set{col}, 'FontSize', 18);
        end
        if col == 1
            ylabel(ylabel_set{row}, 'FontSize', 18);
        end

        xlim([0 1]);
        ylim([0 1]);
        axis square;
        ax = gca;
        set(ax, 'XTick', [0 0.5 1], 'FontSize', 20);
        set(ax, 'YTick', [0 0.5 1], 'FontSize', 20);
        set(ax, 'TickDir', 'out');
    end
end
hgexport(f, [out_dir '/10_fold_CV.eps']);
eps2xxx([out_dir '/10_fold_CV.eps'], {'png'});

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])