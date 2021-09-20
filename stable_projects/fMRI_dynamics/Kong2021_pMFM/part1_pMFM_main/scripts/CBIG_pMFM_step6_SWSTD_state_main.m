function CBIG_pMFM_step6_SWSTD_state_main()

% This function is used to conduct the statistical analysis corresponding
% to figure 5D and 5E in the paper.
%
% There is no input for this function as it can automatically get the
% output file from previous step.
% There is no output for this function as it will generate the output
% files.
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load('../../input/Desikan_input/run_label_testset.mat', 'run_label')
output_dir = '../output/step6_SWSTD_state';
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end

%% Compute empirical SWSTD state
load('../output/step5_STDFCD_results/STD_FCD_empirical_rundata.mat', 'FCD_emp_allrun', 'SWSTD_emp_allrun')
load('../output/step6_SWSTD_state/FCD_threshold_empirical.mat','th_all_emp')

FCD_std = std(FCD_emp_allrun,1,2);
run_num = size(FCD_std,1);
dim = size(SWSTD_emp_allrun,3);

var_up_all = zeros(dim,run_num);
var_down_all = zeros(dim,run_num);

for i = 1:run_num
    FCD_mean = FCD_emp_allrun(i,:);
    threshold = th_all_emp(i);
   
    if threshold ~= 0
        var_time = squeeze(SWSTD_emp_allrun(i,:,:))';
        var_time_demean = bsxfun(@minus, var_time, mean(var_time,2));
        var_time_norm = bsxfun(@rdivide, var_time_demean, std(var_time,1,2));
        var_up = var_time_norm(:,(FCD_mean>threshold));
        var_down = var_time_norm(:,(FCD_mean<threshold));
        var_up_mean = mean(var_up,2);
        var_down_mean = mean(var_down,2);

        var_up_all(:,i) = var_up_mean;
        var_down_all(:,i) = var_down_mean;
    end
     
    disp(['Finish empirical run: ' num2str(i)])
end

sub_num = max(run_label);
sub_up_all = zeros(dim,sub_num);
sub_down_all = zeros(dim,sub_num);
for j = 1:sub_num
    sub_up = var_up_all(:,(run_label==j));
    sub_up(sub_up(1, :)==0) = [];
    sub_up_all(:,j) = mean(sub_up, 2);
    sub_down = var_down_all(:,(run_label==j));
    sub_down(sub_down(1, :)==0) = [];
    sub_down_all(:,j) = mean(sub_down, 2);
end

SWSTD_state_emp_up = mean(sub_up_all,1);
SWSTD_state_emp_down = mean(sub_down_all,1);
SWSTD_state_emp = [SWSTD_state_emp_up', SWSTD_state_emp_down'];
save(fullfile(output_dir, 'SWSTD_state_empirical.mat'), 'SWSTD_state_emp')


%% Compute simulated SWSTD state
load('../output/step5_STDFCD_results/STD_FCD_simulated_rundata.mat', 'FCD_sim_allrun', 'SWSTD_sim_allrun')
load('../output/step6_SWSTD_state/FCD_threshold_simulated.mat','th_all_sim')

FCD_std = std(FCD_sim_allrun,1,2);
run_num = size(FCD_std,1);
dim = size(SWSTD_sim_allrun,3);

var_up_all = zeros(dim,run_num);
var_down_all = zeros(dim,run_num);

for i = 1:run_num
    FCD_mean = FCD_sim_allrun(i,:);
    threshold = th_all_sim(i);
    
    var_time = squeeze(SWSTD_sim_allrun(i,:,:))';
    var_time_demean = bsxfun(@minus, var_time, mean(var_time,2));
    var_time_norm = bsxfun(@rdivide, var_time_demean, std(var_time,1,2));
    var_up = var_time_norm(:,(FCD_mean>threshold));
    var_down = var_time_norm(:,(FCD_mean<threshold));
    var_up_mean = mean(var_up,2);
    var_down_mean = mean(var_down,2);
    
    var_up_all(:,i) = var_up_mean;
    var_down_all(:,i) = var_down_mean;
    
    disp(['Finish simulated run: ' num2str(i)])
end

SWSTD_state_sim_up = mean(var_up_all,1);
SWSTD_state_sim_down = mean(var_down_all,1);
SWSTD_state_sim = [SWSTD_state_sim_up', SWSTD_state_sim_down'];
save(fullfile(output_dir, 'SWSTD_state_simulated.mat'), 'SWSTD_state_sim')

end