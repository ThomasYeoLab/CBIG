function CBIG_pMFM_STDFCD_lowGS()

% This function is used to generate the SWSTD FCD correlation with 50 runs
% with lowest gloable signal
%
% There is no output for this function as it will generate the output files
% which will be used as input for next step
%
% Args:
%    emp_tc_dir: Directory contains empirical fMRI signal files. Data in
%                each file should be in the formate of NxT, where N is 
%                number of ROIs T is the number of time points
%    emp_fcd_dir:Directory contains empirical FCD files. The order of FCD
%                files should be the same as the order of fMRI signal
%                files.
%
% Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
emp_tc_dir = '../../../part0_pMFM_data_preparation/Desikan/output/TC/';
emp_fcd_dir = '../../../part0_pMFM_data_preparation/Desikan/output/FCD/';

stdfcd_lowGS_dir = fullfile('../output');
if ~exist(stdfcd_lowGS_dir,'dir')
    mkdir(stdfcd_lowGS_dir)
end

load('../../../input/Desikan_input/GS_sort.mat', 'data_sort')


window_len = 83;
swstd_all = zeros(68,50);

for i = 1:50
    
    if data_sort(i,3) == 1
        file_index = [num2str(data_sort(i,1)) '_rfMRI_REST' num2str(data_sort(i,2)) '_LR.mat'];
    else
        file_index = [num2str(data_sort(i,1)) '_rfMRI_REST' num2str(data_sort(i,2)) '_RL.mat'];
    end
    
    if isfile(fullfile(emp_tc_dir, 'train', file_index))
        load(fullfile(emp_tc_dir, 'train', file_index))
        load(fullfile(emp_fcd_dir, 'train', file_index))
    elseif isfile(fullfile(emp_tc_dir, 'validation', file_index))
        load(fullfile(emp_tc_dir, 'validation', file_index))
        load(fullfile(emp_fcd_dir, 'validation', file_index))
    elseif isfile(fullfile(emp_tc_dir, 'test', file_index))
        load(fullfile(emp_tc_dir, 'test', file_index))
        load(fullfile(emp_fcd_dir, 'test', file_index))
    else
        continue
    end
        
    if size(FCD_mat,1) < 1118
        continue;
    end
    
    TC([1,5,37,41],:) = [];
    FCD_mean = mean(FCD_mat,1);
    FCD_grad = FCD_mean(:,2:end)-FCD_mean(:,1:end-1);
    
    var_time = zeros(68,1200-window_len+1,1);
    for ti = 1:length(var_time)
        var_time(:,ti) = std(TC(:,ti:ti+window_len-1),1,2);
    end
    
    var_grad = var_time(:,2:end)-var_time(:,1:end-1);
    
    corr_all = corr(var_grad',FCD_grad','type','Spearman');
    swstd_all(:,i) = corr_all;

end
swstd_all_nz = swstd_all;
swstd_all_nz(:,swstd_all_nz(1,:)==0) = [];
SWSTD_FCD_emp = nanmean(swstd_all_nz,2);
save(fullfile(stdfcd_lowGS_dir, 'STD_FCD_lowGS_empirical.mat'), 'SWSTD_FCD_emp')
    
end

    
