function CBIG_pMFM_step1_generate_permutation_order_desikan()

% This function is the first step of the permutation test of SWSTD-FCD
% correlation. In this step, 10000 permutation order will be generated and
% store for further use.
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% generate permutation order for empirical data
disp('Start generating empirical order')
for e = 1:100
    tic
    generate_permutation_order_empirical(e)
    toc
    disp(e)
end

%% generate permutation order for simulated data
disp('Start generating simulated order')
for s = 1:100
    tic
    generate_permutation_order_simulated(s)
    toc
    disp(s)
end
end


function generate_permutation_order_empirical(o)

% This function is used to generate permutation order for empirical data
% In total we perform 10000 permuations
% Args:
%    o: random seed for 100 permute
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load('../input/run_label_testset.mat','run_label');
run_num = size(run_label, 1);
perm_time = 100;
perm_result_all = zeros(68,run_num,perm_time);
rng(o)
for r = 1:68
   
    perm_result = zeros(run_num,perm_time);
    p = 1;
    while p <= perm_time
        perm_label = 1:run_num;
        perm_indi = ones(run_num,1);
        for i = 1:run_num
            temp_indi = (run_label~=run_label(i))*1.*perm_indi;
            if sum(temp_indi) ~= 0
                temp_label = perm_label(temp_indi==1);
                rand_num = randi(length(temp_label),1);
                perm_result(i,p) = temp_label(rand_num);
                perm_indi(temp_label(rand_num)) = 0;
            else
                p = p-1;
                break
            end
        end
        p = p+1;
    end
    perm_result_all(r,:,:) = perm_result;
end

emp_dir = '../output/step1_generate_order/empirical';
if ~exist(emp_dir, 'dir')
    mkdir(emp_dir)
end

save([emp_dir '/perm_order_' num2str(o) '.mat'],'perm_result_all')

end


function generate_permutation_order_simulated(o)

% This function is used to generate permutation order for simulated data
% In total we perform 10000 permuations
% Args:
%    o: random seed for 100 permute
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load('../../../part1_pMFM_main/output/step5_STDFCD_results/STD_FCD_simulated_rundata.mat','FCD_sim_allrun')
run_num = size(FCD_sim_allrun,1);

rng(o)
perm_time = 100;
perm_result_all = zeros(68,run_num, perm_time);


for r = 1:68
    perm_result = zeros(run_num,perm_time);
    for i = 1:perm_time
        perm_result(:, i) = randperm(run_num);
    end
    perm_result_all(r,:,:) = perm_result;
end

sim_dir = '../output/step1_generate_order/simulated';
if ~exist(sim_dir, 'dir')
    mkdir(sim_dir)
end

save([sim_dir '/perm_order_' num2str(o) '.mat'], 'perm_result_all')


end


