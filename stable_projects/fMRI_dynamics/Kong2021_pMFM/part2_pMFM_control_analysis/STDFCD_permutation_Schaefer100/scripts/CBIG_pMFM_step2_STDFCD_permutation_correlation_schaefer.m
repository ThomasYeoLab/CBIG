function CBIG_pMFM_step2_STDFCD_permutation_correlation()

% This function is the second step of the permutation test of SWSTD-FCD
% correlation. In this step, 10000 permuted SWSTD FCD correlation will be 
% generated and store for further use.
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% generate permuted correlation for empirical data
disp('Start generating empirical correlation')
for e = 1:100
    tic
    compute_STDFCD_correlation_empirical(e)
    toc
    disp(e)
end

%% generate permutation order for simulated data
disp('Start generating simulated correlation')
for s = 1:100
    tic
    compute_STDFCD_correlation_simulated(s)
    toc
    disp(s)
end

%% Perform statistical test
statistical_test_empirical()
statistical_test_simulated()

end



function compute_STDFCD_correlation_empirical(o)

% This function loads the permutation order and compute SWSTD-FCD
% correlation to generate null distribution for empirical data.
% Args:
%    o: random seed for 100 permute
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load('../../Schaefer100_parcellation/output/step5_STDFCD_results/STD_FCD_empirical_rundata.mat', ...
    'FCD_emp_allrun', 'SWSTD_emp_allrun');

load(['../output/step1_generate_order/empirical/perm_order_' num2str(o) '.mat'],'perm_result_all');

sub_num = size(FCD_emp_allrun,1);
perm_time = size(perm_result_all,3);
null_mat = zeros(100,perm_time);

for i = 1:perm_time
    perm_index = perm_result_all(:,:,i);
    corr_perm = zeros(100,sub_num);
    for j = 1:100
        var_grad = SWSTD_emp_allrun(:,:,j);
        fcd_grad_perm = FCD_emp_allrun(perm_index(j,:),:);
        corr_ROI = diag(corr(var_grad',fcd_grad_perm','type','Spearman'));
        corr_perm(j,:) = corr_ROI;
    end
    null_mat(:,i) = mean(corr_perm,2);
end
emp_dir = '../output/step2_compute_correlation/empirical';
if ~exist(emp_dir, 'dir')
    mkdir(emp_dir)
end

save([emp_dir '/correlation_' num2str(o) '.mat'],'null_mat')

end



function compute_STDFCD_correlation_simulated(o)

% This function loads the permutation order and compute SWSTD-FCD
% correlation to generate null distribution for simulated data.
% Args:
%    o: random seed for 100 permute
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load('../../Schaefer100_parcellation/output/step5_STDFCD_results/STD_FCD_simulated_rundata.mat', ...
    'FCD_sim_allrun','SWSTD_sim_allrun');
load(['../output/step1_generate_order/simulated/perm_order_' num2str(o) '.mat'],'perm_result_all');

sub_num = size(FCD_sim_allrun,1);
perm_time = size(perm_result_all,3);
null_mat = zeros(100,perm_time);

for i = 1:perm_time
    perm_index = perm_result_all(:,:,i);
    corr_perm = zeros(100,sub_num);
    for j = 1:100
        var_grad = SWSTD_sim_allrun(:,:,j);
        fcd_grad_perm = FCD_sim_allrun(perm_index(j,:),:);
        corr_ROI = diag(corr(var_grad',fcd_grad_perm','type','Spearman'));
        corr_perm(j,:) = corr_ROI;
    end
    null_mat(:,i) = mean(corr_perm,2);
end

sim_dir = '../output/step2_compute_correlation/simulated';
if ~exist(sim_dir, 'dir')
    mkdir(sim_dir)
end

save([sim_dir '/correlation_' num2str(o) '.mat'],'null_mat')

end


function statistical_test_empirical()

% This function loads the computed SWSTD-FCD correlations and perform
% statistical test.
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

emp_dir = '../output/step2_compute_correlation/empirical';
null_test_file = dir(emp_dir);
null_test_all = zeros(100,100*(length(null_test_file)-2));

load('../../Schaefer100_parcellation/output/step5_STDFCD_results/STD_FCD_empirical.mat','SWSTD_FCD_emp')

for i = 1:length(null_test_file)-2
    load(fullfile(emp_dir, null_test_file(i+2).name),'null_mat')
    null_test_all(:,500*(i-1)+1:500*i) = null_mat;
end

p_value = zeros(100,1);
null_abs = abs(null_test_all);
for i = 1:100
    p_value(i) = (sum(null_abs(i,:)>SWSTD_FCD_emp(i))+1)/10001;
end

save('output/STDFCD_correlation_pvalue_empirical.mat','p_value')

end


function statistical_test_simulated()

% This function loads the computed SWSTD-FCD correlations and perform
% statistical test.
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

sim_dir = '../output/step2_compute_correlation/simulated';
null_test_file = dir(sim_dir);
null_test_all = zeros(100,100*(length(null_test_file)-2));

load('../../Schaefer100_parcellation/output/step5_STDFCD_results/STD_FCD_simulated.mat','SWSTD_FCD_sim')

for i = 1:length(null_test_file)-2
    load(fullfile(sim_dir, null_test_file(i+2).name),'null_mat')
    null_test_all(:,500*(i-1)+1:500*i) = null_mat;
end

p_value = zeros(100,1);
null_abs = abs(null_test_all);
for i = 1:100
    p_value(i) = (sum(null_abs(i,:)>SWSTD_FCD_sim(i))+1)/10001;
end

save('output/STDFCD_correlation_pvalue_simulated.mat','p_value')

end
        
