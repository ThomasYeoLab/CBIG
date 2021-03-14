function CBIG_pMFM_step5_generate_STDFCD_correlation_main()

% This function is used to generate the correlation between the empirical
% and simualted SWSTD of timecourse. The results corresponds to the figure
% 6A, 6B and 6C in the paper.
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output files
% which will be used as input for next step
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

stdfcd_result_dir = fullfile('../output/step5_STDFCD_results');
if ~exist(stdfcd_result_dir,'dir')
    mkdir(stdfcd_result_dir)
end


%% Empirical SWSTD generation
emp_tc_dir = '../../part0_pMFM_data_preparation/Desikan/output/TC/test';
emp_fcd_dir = '../../part0_pMFM_data_preparation/Desikan/output/FCD/test';

emp_tc_file = dir(emp_tc_dir);
emp_fcd_file = dir(emp_fcd_dir);

window_len = 83;
file_num = length(emp_fcd_file)-2;
swstd_all = zeros(68,file_num);

FCD_emp_allrun = zeros(file_num,1200-window_len);
SWSTD_emp_allrun = zeros(file_num,1200-window_len,68);

for i = 3:file_num+2
    load(fullfile(emp_tc_dir, emp_tc_file(i).name), 'TC')
    load(fullfile(emp_fcd_dir, emp_fcd_file(i).name), 'FCD_mat')

    FCD_length = size(FCD_mat,1);
    if FCD_length < 1200-window_len+1
        continue;
    end

    TC([1,5,37,41],:) = [];
    FCD_mean = mean(FCD_mat,1);
    FCD_grad = FCD_mean(:,2:end)-FCD_mean(:,1:end-1);
    FCD_emp_allrun(i-2, :) = FCD_grad;

    var_time = zeros(68,1200-window_len+1,1);
    for ti = 1:length(var_time)
        var_time(:,ti) = std(TC(:,ti:ti+window_len-1),1,2);
    end
    var_grad = var_time(:,2:end)-var_time(:,1:end-1);
    SWSTD_emp_allrun(i-2,:,:) = var_grad';

    TC_grad_sum = corr(var_grad',FCD_grad','type','Spearman');

    swstd_all(:,i-2) = TC_grad_sum;
    disp(['Finish subject: ' num2str(i-2)])

end

swstd_all_nz = swstd_all;
swstd_all_nz(:,swstd_all_nz(1,:)==0) = [];
SWSTD_FCD_emp = nanmean(swstd_all_nz,2);
SWSTD_FCD_emp_all = swstd_all_nz;
% Save out averaged empirical SWSTD-FCD correlation into SWSTD_FCD_emp
% Save out all the empirical SWSTD-FCD correlations into SWSTD_FCD_emp_all
save(fullfile(stdfcd_result_dir, 'STD_FCD_empirical.mat'), 'SWSTD_FCD_emp','SWSTD_FCD_emp_all')

FCD_emp_allrun_sum = sum(FCD_emp_allrun,2);
FCD_emp_allrun(FCD_emp_allrun_sum==0,:) = [];
SWSTD_emp_allrun_sum = sum(sum(SWSTD_emp_allrun,3),2);
SWSTD_emp_allrun(SWSTD_emp_allrun_sum==0,:,:) = [];
% Save out all the empirical FCD mean and SWSTD for next step use
save(fullfile(stdfcd_result_dir, 'STD_FCD_empirical_rundata.mat'), 'FCD_emp_allrun', 'SWSTD_emp_allrun')


%% Simulated SWSTD generation
sim_tc_dir = '../output/step4_MFM_simulated_data/TC';
sim_tc_file = dir(sim_tc_dir);
sim_fcd_dir = '../output/step4_MFM_simulated_data/FCD';
sim_fcd_file = dir(sim_fcd_dir);

window_len = 83;
file_num = length(sim_tc_file)-2;
swstd_all = zeros(68,file_num);

FCD_sim_allrun = zeros(file_num,1200-window_len);
SWSTD_sim_allrun = zeros(file_num,1200-window_len,68);

for i = 3:file_num+2
    load(fullfile(sim_tc_dir, sim_tc_file(i).name), 'TC')
    load(fullfile(sim_fcd_dir, sim_fcd_file(i).name), 'FCD_mat')

    FCD_length = size(FCD_mat,1);
    if FCD_length < 1200-window_len+1
        continue;
    end

    FCD_mean = mean(FCD_mat,1);
    FCD_grad = FCD_mean(:,2:end)-FCD_mean(:,1:end-1);
    FCD_sim_allrun(i-2, :) = FCD_grad;

    var_time = zeros(68,1200-window_len+1,1);
    for ti = 1:length(var_time)
        var_time(:,ti) = std(TC(:,ti:ti+window_len-1),1,2);
    end
    var_grad = var_time(:,2:end)-var_time(:,1:end-1);
    SWSTD_sim_allrun(i-2,:,:) = var_grad';

    TC_grad_sum = corr(var_grad',FCD_grad','type','Spearman');

    swstd_all(:,i-2) = TC_grad_sum;
    disp(['Finish subject: ' num2str(i-2)])

end

swstd_all_nz = swstd_all;
swstd_all_nz(:,swstd_all_nz(1,:)==0) = [];
SWSTD_FCD_sim = nanmean(swstd_all_nz,2);
SWSTD_FCD_sim_all = swstd_all_nz;
% Save out averaged simulated SWSTD-FCD correlation into SWSTD_FCD_sim
% Save out all the simulated SWSTD-FCD correlations into SWSTD_FCD_sim_all
save(fullfile(stdfcd_result_dir, 'STD_FCD_simulated.mat'), 'SWSTD_FCD_sim','SWSTD_FCD_sim_all')

FCD_sim_allrun_sum = sum(FCD_sim_allrun,2);
FCD_sim_allrun(FCD_sim_allrun_sum==0,:) = [];
SWSTD_sim_allrun_sum = sum(sum(SWSTD_sim_allrun,3),2);
SWSTD_sim_allrun(SWSTD_sim_allrun_sum==0,:,:) = [];
% Save out all the simulated FCD mean and SWSTD for next step use
save(fullfile(stdfcd_result_dir, 'STD_FCD_simulated_rundata.mat'), 'FCD_sim_allrun', 'SWSTD_sim_allrun')

%% Compute correlation
STDFCD_emp_sim_corr = corr(SWSTD_FCD_emp, SWSTD_FCD_sim,'type','Spearman');
save(fullfile(stdfcd_result_dir, 'STD_FCD_correlation.mat'), 'STDFCD_emp_sim_corr')

end