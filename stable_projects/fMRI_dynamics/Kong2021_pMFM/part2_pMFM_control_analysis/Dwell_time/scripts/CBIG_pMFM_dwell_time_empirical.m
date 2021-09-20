function CBIG_pMFM_dwell_time_empirical()

% This function is used to compute the dwell time distirbution for empirial
% data
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load('../../../part1_pMFM_main/output/step5_STDFCD_results/STD_FCD_empirical_rundata.mat', ...
    'FCD_emp_allrun')
load('../../../part1_pMFM_main/output/step6_SWSTD_state/FCD_threshold_empirical.mat', ...
    'th_all_emp')

FCD_std = std(FCD_emp_allrun,1,2);
run_num = size(FCD_std,1);

count_up = [];
count_down = [];


for i = 1:run_num
    FCD_mean = FCD_emp_allrun(i,:);
    threshold = th_all_emp(i);
    
    if threshold ~= 0
        FCD_up = 1*(FCD_mean>threshold);
        FCD_down = 1*(FCD_mean<threshold);

        FCD_up_count = CBIG_pMFM_count_func(FCD_up);
        FCD_down_count = CBIG_pMFM_count_func(FCD_down);

        count_up = [count_up;FCD_up_count];
        count_down = [count_down;FCD_down_count];
    end
end

[dwell_up_emp, edges_up] = histcounts(count_up, 1:12:1201);
bar(log(dwell_up_emp))
[dwell_down_emp, edges_down] = histcounts(count_down, 1:12:1201);
bar(log(dwell_down_emp))

output_dir = '../output';
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end

save('output/dwell_high_emp.mat','dwell_up_emp')
save('output/dwell_down_emp.mat','dwell_down_emp')

end
