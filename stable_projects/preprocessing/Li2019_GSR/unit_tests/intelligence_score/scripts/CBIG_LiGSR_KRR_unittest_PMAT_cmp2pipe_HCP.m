function CBIG_LiGSR_KRR_unittest_PMAT_cmp2pipe_HCP( KRR_dir )

% CBIG_LiGSR_KRR_unittest_PMAT_cmp2pipe_HCP( KRR_dir )
% 
% This is step 2 of the unit test which predicts fluid intelligence score
% in the HCP dataset.
% 
% This function integrates the output files of
% CBIG_LiGSR_KRR_uniittest_PMAT_HCP.sh and compares the results between 
% Baseline and Baseline_GSR pipelines.
% 
% Inputs:
%   - KRR_dir
%     The output directory of CBIG_LiGSR_KRR_uniittest_PMAT_HCP.sh. It is
%     assumed that there are two subfolders call "Baseline" and "GSR" 
%     storing the kernel ridge regression results generated using each
%     preprocessing pipeline.
%     A file called [KRR_dir '/compare_2pipe/final_result.mat'] will be
%     created to store the integrated result and the comparison between the
%     two pipelines. It contains:
%     (1) mean_acc_dif: a 20x1 vector of fluid intelligence prediction
%                       accuracy difference between the two preprocessing
%                       pipelines, averaged within each random data split.
%     (2) acc_GSR_mean: a 20x1 vector of accuracy of fluid intelligence
%                       with global signal regression, averaged within each
%                       random data split.
%     (3) acc_Baseline_mean: a 20x1 vector of accuracy of fluid
%                            intelligence with baseline fMRI preprocessing,
%                            averaged within each random data split.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~exist(fullfile(KRR_dir, 'compare_2pipe'), 'dir'))
    mkdir(fullfile(KRR_dir, 'compare_2pipe'))
end

stem = {'PMAT24'};

for seed = 1:3
    GSR = load(fullfile(KRR_dir, 'GSR', ['randseed_' num2str(seed)], ...
        ['final_result_' stem{1} '.mat']));
    Baseline = load(fullfile(KRR_dir, 'Baseline', ['randseed_' num2str(seed)], ...
        ['final_result_' stem{1} '.mat']));
        
    acc_GSR_mean(seed, :) = mean(GSR.optimal_acc, 1);
    acc_Baseline_mean(seed, :) = mean(Baseline.optimal_acc, 1);
end

mean_acc_dif = acc_GSR_mean - acc_Baseline_mean;
save(fullfile(KRR_dir, 'compare_2pipe', 'final_result.mat'), 'mean_acc_dif', 'acc_GSR_mean', ...
    'acc_Baseline_mean')

end
