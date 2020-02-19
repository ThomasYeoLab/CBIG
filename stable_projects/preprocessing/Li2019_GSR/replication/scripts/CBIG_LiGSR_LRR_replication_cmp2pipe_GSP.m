function CBIG_LiGSR_LRR_replication_cmp2pipe_GSP( LRR_dir )

% CBIG_LiGSR_LRR_replication_cmp2pipe_GSP( LRR_dir )
% 
% This is step 2 of the replication unit test using the GSP dataset.
% 
% This function integrate the output files of
% CBIG_LiGSR_LRR_replication_all_GSP.sh and compare the results between 
% Baseline and Baseline_GSR pipelines.
% 
% Inputs:
%   - LRR_dir
%     The output directory of CBIG_LiGSR_LRR_replication_all_GSP.sh. It is
%     assumed that there are two subfolders call "Baseline" and "GSR" storing
%     the linear ridge regression results generated using each
%     preprocessing pipeline.
% 
% Outputs:
%     A file called [LRR_dir '/compare_2pipe/final_result.mat'] will be
%     created to store the integrated result and the comparison between the
%     two pipelines. It contains:
%     (1) mean_acc_dif: a #seed x #MeasuresToPredict (i.e. 20 x 24) matrix
%                       of accuracy difference between the two
%                       preprocessing pipelines, averaged within each
%                       random data split.
%     (2) acc_GSR_mean: a #seed x #MeasuresToPredict (i.e. 20 x 24) matrix
%                       of accuracy with global signal regression, averaged
%                       within each random data split.
%     (3) acc_Baseline_mean: a #seed x #MeasuresToPredict (i.e. 20 x 24)
%                            matrix of accuracy with baseline fMRI
%                            preprocessing, averaged within each random
%                            data split.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~exist(fullfile(LRR_dir, 'compare_2pipe'), 'dir'))
    mkdir(fullfile(LRR_dir, 'compare_2pipe'))
end

behavior_txt = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'replication', 'scripts', 'GSP_lists', '23behaviors.txt');
Age_txt = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', 'Li2019_GSR', ...
    'replication', 'scripts', 'GSP_lists', 'Age_header.txt');
behavior_stem = CBIG_text2cell(behavior_txt);
Age_stem = CBIG_text2cell(Age_txt);

stem = [behavior_stem Age_stem];

for seed = 1:20
    for i = 1:numel(stem)
        GSR = load(fullfile(LRR_dir, 'GSR', ['randseed_' num2str(seed)], ...
            'results', 'optimal_acc', [stem{i} '.mat']));
        Baseline = load(fullfile(LRR_dir, 'Baseline', ['randseed_' num2str(seed)], ...
            'results', 'optimal_acc', [stem{i} '.mat']));
        
        if(i == 1)
            opt_GSR_mean_tmp = mean(GSR.acc_corr_test, 1);
            opt_Base_mean_tmp = mean(Baseline.acc_corr_test, 1);
        else            
            opt_GSR_mean_tmp = [opt_GSR_mean_tmp mean(GSR.acc_corr_test, 1)];
            opt_Base_mean_tmp = [opt_Base_mean_tmp mean(Baseline.acc_corr_test, 1)];
        end
    end
    acc_GSR_mean(seed, :) = opt_GSR_mean_tmp;
    acc_Baseline_mean(seed, :) = opt_Base_mean_tmp;
end

mean_acc_dif = acc_GSR_mean - acc_Baseline_mean;
save(fullfile(LRR_dir, 'compare_2pipe', 'final_result.mat'), 'mean_acc_dif', ...
    'acc_GSR_mean', 'acc_Baseline_mean')

end

