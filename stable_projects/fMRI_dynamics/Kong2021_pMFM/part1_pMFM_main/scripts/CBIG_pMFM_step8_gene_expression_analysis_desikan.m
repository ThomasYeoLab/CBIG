function CBIG_pMFM_step8_gene_expression_analysis_desikan()

% This function is the wrapper to perform all the genetic analysis
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

gene_result_dir = '../output/step8_gene_expression_analysis';
if ~exist(gene_result_dir,'dir')
    mkdir(gene_result_dir)
end

gene_MFMresult_correlation()
gene_bootstrapping_analysis()

end

function gene_MFMresult_correlation()

% This function is used to generate the correlation between the gene
% expression data (e.g. PVALB-SST) and MFM analysis results (e.g estmiated 
% model parameters) and compute the p-value.
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
% which will be used as input for next step
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% MFM estimated model parameters
para_all = csvread('../output/step3_test_results/test_all.csv');
para_first = para_all(12:end, 1);
w_para = para_first(1:68);
I_para = para_first(69:136);
s_para = para_first(138:end);


%% SWSTD-FCD correlation data
load('../output/step5_STDFCD_results/STD_FCD_empirical.mat', 'SWSTD_FCD_emp')
load('../output/step5_STDFCD_results/STD_FCD_simulated.mat', 'SWSTD_FCD_sim')
    
%% gene expression data
load('../../input/Desikan_input/gene_data.mat','PVALB','SST','PC1')
gene_diff = PVALB-SST;

%% compute correlation
PVALBSST_STDFCD_correlation_empirical = corr(SWSTD_FCD_emp, gene_diff, 'type', 'Spearman');
PC1_STDFCD_correlation_empirical = corr(SWSTD_FCD_emp, PC1, 'type', 'Spearman');
PVALBSST_STDFCD_correlation_simulated = corr(SWSTD_FCD_sim, gene_diff, 'type', 'Spearman');
PC1_STDFCD_correlation_simulated = corr(SWSTD_FCD_sim, PC1, 'type', 'Spearman');

PVALBSST_w_correlation = corr(w_para, gene_diff, 'type', 'Spearman');
PC1_w_correlation = corr(w_para, PC1, 'type', 'Spearman');
PVALBSST_I_correlation = corr(I_para, gene_diff, 'type', 'Spearman');
PC1_I_correlation = corr(I_para, PC1, 'type', 'Spearman');
PVALBSST_s_correlation = corr(s_para, gene_diff, 'type', 'Spearman');
PC1_s_correlation = corr(s_para, PC1, 'type', 'Spearman');


save('../output/step8_gene_expression_analysis/gene_expression_MFM_correlation.mat',...
    'PVALBSST_STDFCD_correlation_empirical','PC1_STDFCD_correlation_empirical',...
    'PVALBSST_STDFCD_correlation_simulated','PC1_STDFCD_correlation_simulated',...
    'PVALBSST_w_correlation','PC1_w_correlation','PVALBSST_I_correlation',...
    'PC1_I_correlation','PVALBSST_s_correlation','PC1_s_correlation')


%% load spin test data
load('../../input/Desikan_input/gene_spin_data.mat','PC1_spin_data','PVALB_spin_data','SST_spin_data')
gene_diff_spin = PVALB_spin_data-SST_spin_data;

%% compute p-value for spin test
PVALBSST_STDFCD_spin_empirical = corr(gene_diff_spin', SWSTD_FCD_emp, 'type', 'Spearman');
p_PVALBSST_STDFCD_empirical = (sum(abs(PVALBSST_STDFCD_spin_empirical)>...
    abs(PVALBSST_STDFCD_correlation_empirical))+1)/1001;

PVALBSST_STDFCD_spin_simulated = corr(gene_diff_spin', SWSTD_FCD_sim, 'type', 'Spearman');
p_PVALBSST_STDFCD_simulated = (sum(abs(PVALBSST_STDFCD_spin_simulated)>...
    abs(PVALBSST_STDFCD_correlation_simulated))+1)/1001;

PC1_STDFCD_spin_empirical = corr(PC1_spin_data', SWSTD_FCD_emp, 'type', 'Spearman');
p_PC1_STDFCD_empirical = (sum(abs(PC1_STDFCD_spin_empirical)>...
    abs(PC1_STDFCD_correlation_empirical))+1)/1001;

PC1_STDFCD_spin_simulated = corr(PC1_spin_data', SWSTD_FCD_sim, 'type', 'Spearman');
p_PC1_STDFCD_simulated = (sum(abs(PC1_STDFCD_spin_simulated)>...
    abs(PC1_STDFCD_correlation_simulated))+1)/1001;

PVALBSST_w_spin = corr(gene_diff_spin', w_para, 'type', 'Spearman');
p_PVALBSST_w = (sum(abs(PVALBSST_w_spin)>abs(PVALBSST_w_correlation))+1)/1001;

PVALBSST_I_spin = corr(gene_diff_spin', I_para, 'type', 'Spearman');
p_PVALBSST_I = (sum(abs(PVALBSST_I_spin)>abs(PVALBSST_I_correlation))+1)/1001;

PVALBSST_s_spin = corr(gene_diff_spin', s_para, 'type', 'Spearman');
p_PVALBSST_s = (sum(abs(PVALBSST_s_spin)>abs(PVALBSST_s_correlation))+1)/1001;

PC1_w_spin = corr(PC1_spin_data', w_para, 'type', 'Spearman');
p_PC1_w = (sum(abs(PC1_w_spin)>abs(PC1_w_correlation))+1)/1001;

PC1_I_spin = corr(PC1_spin_data', I_para, 'type', 'Spearman');
p_PC1_I = (sum(abs(PC1_I_spin)>abs(PC1_I_correlation))+1)/1001;

PC1_s_spin = corr(PC1_spin_data', s_para, 'type', 'Spearman');
p_PC1_s = (sum(abs(PC1_s_spin)>abs(PC1_s_correlation))+1)/1001;

save('../output/step8_gene_expression_analysis/gene_expression_MFM_pvalue_spin_test.mat',...
    'p_PVALBSST_STDFCD_empirical','p_PVALBSST_STDFCD_simulated',...
    'p_PC1_STDFCD_empirical','p_PC1_STDFCD_simulated',...
    'p_PVALBSST_w','p_PVALBSST_I','p_PVALBSST_s',...
    'p_PC1_w','p_PC1_I','p_PC1_s')

%% load random gene data
gene_random_data = csvread('../../input/Desikan_input/gene_expression_data.csv');
gene_random_data(:,[4,39]) = [];

%% compute p-value for random gene test
gene_num = size(gene_random_data, 1);
rng(0)
random_num = randi([1, gene_num], 10000, 2);
gene_random_diff = zeros(10000, size(gene_random_data, 2));
for i = 1:10000
    gene_random_diff(i, :) = gene_random_data(random_num(i, 1), :)-gene_random_data(random_num(i, 2), :);
end


PVALBSST_STDFCD_random_empirical = corr(gene_random_diff', SWSTD_FCD_emp, 'type', 'Spearman');
p_PVALBSST_STDFCD_empirical_random = (sum(abs(PVALBSST_STDFCD_random_empirical)>...
    abs(PVALBSST_STDFCD_correlation_empirical))+1)/10001;

PVALBSST_STDFCD_random_simulated = corr(gene_random_diff', SWSTD_FCD_sim, 'type', 'Spearman');
p_PVALBSST_STDFCD_simulated_random = (sum(abs(PVALBSST_STDFCD_random_simulated)>...
    abs(PVALBSST_STDFCD_correlation_simulated))+1)/10001;

PVALBSST_w_random = corr(gene_random_diff', w_para, 'type', 'Spearman');
p_PVALBSST_w_random = (sum(abs(PVALBSST_w_random)>abs(PVALBSST_w_correlation))+1)/10001;

PVALBSST_I_random = corr(gene_random_diff', I_para, 'type', 'Spearman');
p_PVALBSST_I_random = (sum(abs(PVALBSST_I_random)>abs(PVALBSST_I_correlation))+1)/10001;

PVALBSST_s_random = corr(gene_random_diff', s_para, 'type', 'Spearman');
p_PVALBSST_s_random = (sum(abs(PVALBSST_s_random)>abs(PVALBSST_s_correlation))+1)/10001;

save('../output/step8_gene_expression_analysis/gene_expression_MFM_pvalue_random_test.mat',...
    'p_PVALBSST_STDFCD_empirical_random','p_PVALBSST_STDFCD_simulated_random',...
    'p_PVALBSST_w_random','p_PVALBSST_I_random','p_PVALBSST_s_random')


end

function gene_bootstrapping_analysis()

% This function is used to generate the SWSTD-FCD correlation distribution
% by using boostrapping method and then correlate with gene expression data
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% bootstrapping on empirical SWSTD-FCD
load('../../input/Desikan_input/run_label_testset.mat','run_label')
rng(1)
sub_num = max(run_label);
run_num = zeros(sub_num,1);
for i = 1:sub_num
    run_num(i) = sum(run_label==i);
end
num_4run = sum(run_num==4);
num_3run = sum(run_num==3);
num_2run = sum(run_num==2);
sub_4run = find(run_num==4);
sub_3run = find(run_num==3);
sub_2run = find(run_num==2);

load('../output/step5_STDFCD_results/STD_FCD_empirical.mat', 'SWSTD_FCD_emp_all')
BT_num = 1000;
BT_correlation_empirical = zeros(BT_num, 68);
for bt = 1:BT_num
    BT_corr_Emp = [];
        for i = 1:num_4run
            sub_sel = sub_4run(randi(num_4run,1));
            corr_sel = SWSTD_FCD_emp_all(:,run_label==sub_sel);
            BT_corr_Emp = [BT_corr_Emp,corr_sel];
        end

        for i = 1:num_3run
            sub_sel = sub_3run(randi(num_3run,1));
            corr_sel = SWSTD_FCD_emp_all(:,run_label==sub_sel);
            BT_corr_Emp = [BT_corr_Emp,corr_sel];
        end

        for i = 1:num_2run
            sub_sel = sub_2run(randi(num_2run,1));
            corr_sel = SWSTD_FCD_emp_all(:,run_label==sub_sel);
            BT_corr_Emp = [BT_corr_Emp,corr_sel];
        end
    BT_mean = mean(BT_corr_Emp,2);
    BT_correlation_empirical(bt, :) = BT_mean';
end


%% bootstrapping on simulated SWSTD-FCD
load('../output/step5_STDFCD_results/STD_FCD_simulated.mat', 'SWSTD_FCD_sim_all')

sim_num = size(SWSTD_FCD_sim_all,2);
BT_num = 1000;
BT_correlation_simulated = zeros(BT_num, 68);
rng(1)

for bt = 1:BT_num
    BT_corr_sim = SWSTD_FCD_sim_all(:,randi(sim_num,sim_num,1));
    BT_mean = mean(BT_corr_sim,2);
    BT_correlation_simulated(bt, :) = BT_mean';
end

%% compute bootstrapping SWSTD-FCD and MFM results correlation
load('../../input/Desikan_input/gene_data.mat','PC1','PVALB','SST')
gene_diff = PVALB-SST;

PVALBSST_BT_empirical_correlation = corr(BT_correlation_empirical', gene_diff, 'type', 'Spearman');
PVALBSST_BT_simulated_correlation = corr(BT_correlation_simulated', gene_diff, 'type', 'Spearman');
PC1_BT_empirical_correlation = corr(BT_correlation_empirical', PC1, 'type', 'Spearman');
PC1_BT_simulated_correlation = corr(BT_correlation_empirical', PC1, 'type', 'Spearman');

save('../output/step8_gene_expression_analysis/gene_expression_MFM_bootstrapping.mat',...
    'PVALBSST_BT_empirical_correlation','PVALBSST_BT_simulated_correlation',...
    'PC1_BT_empirical_correlation','PC1_BT_simulated_correlation')

end



