function CBIG_LiGSR_KRR_replication_cmp2pipe_GSP( KRR_dir )

% CBIG_LiGSR_KRR_replication_cmp2pipe_GSP( KRR_dir )
% 
% This is step 2 of the replication unit test using the GSP dataset.
% 
% This function integrate the output files of
% CBIG_LiGSR_KRR_replication_all_GSP.sh and compare the results between 
% Baseline and Baseline_GSR pipelines.
% 
% Inputs:
%   - KRR_dir
%     The output directory of CBIG_LiGSR_KRR_replication_all_GSP.sh. It is
%     assumed that there are two subfolders call "Baseline" and "GSR" storing
%     the kernel ridge regression results generated using each
%     preprocessing pipeline.
%     A file called [KRR_dir '/compare_2pipe/final_result.mat'] will be
%     created to store the integrated result and the comparison between the
%     two pipelines. It contains:
%     (1) mean_acc_dif: a #seed x #MeasuresToPredict (i.e. 20 x 60) matrix
%                       of accuracy difference between the two
%                       preprocessing pipelines, averaged within each
%                       random data split.
%     (2) acc_GSR_mean: a #seed x #MeasuresToPredict (i.e. 20 x 60) matrix
%                       of accuracy with global signal regression, averaged
%                       within each random data split.
%     (3) acc_Baseline_mean: a #seed x #MeasuresToPredict (i.e. 20 x 60)
%                            matrix of accuracy with baseline fMRI
%                            preprocessing, averaged within each random
%                            data split.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~exist(fullfile(KRR_dir, 'compare_2pipe'), 'dir'))
    mkdir(fullfile(KRR_dir, 'compare_2pipe'))
end

stem = {'23behaviors', 'Age', 'Sex'};

for seed = 1:20
    for i = 1:numel(stem)
        GSR = load(fullfile(KRR_dir, 'GSR', ['randseed_' num2str(seed)], ...
            ['final_result_' stem{i} '.mat']));
        Baseline = load(fullfile(KRR_dir, 'Baseline', ['randseed_' num2str(seed)], ...
            ['final_result_' stem{i} '.mat']));
        
        if(i == 1)
            opt_GSR_mean_tmp = mean(GSR.optimal_acc, 1);
            opt_Base_mean_tmp = mean(Baseline.optimal_acc, 1);
        else            
            opt_GSR_mean_tmp = [opt_GSR_mean_tmp mean(GSR.optimal_acc, 1)];
            opt_Base_mean_tmp = [opt_Base_mean_tmp mean(Baseline.optimal_acc, 1)];
        end
    end
    acc_GSR_mean(seed, :) = opt_GSR_mean_tmp;
    acc_Baseline_mean(seed, :) = opt_Base_mean_tmp;
end

mean_acc_dif = acc_GSR_mean - acc_Baseline_mean;
save(fullfile(KRR_dir, 'compare_2pipe', 'final_result.mat'), 'mean_acc_dif', ...
    'acc_GSR_mean', 'acc_Baseline_mean')

end

