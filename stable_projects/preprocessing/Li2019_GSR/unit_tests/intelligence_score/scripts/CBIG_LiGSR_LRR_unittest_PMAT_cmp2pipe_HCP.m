function CBIG_LiGSR_LRR_unittest_PMAT_cmp2pipe_HCP( LRR_dir )

% CBIG_LiGSR_LRR_unittest_PMAT_cmp2pipe_HCP( LRR_dir )
% 
% This is step 2 of the unit test which predicts fluid intelligence score
% in the HCP dataset.
% 
% This function integrates the output files of
% CBIG_LiGSR_LRR_uniittest_PMAT_HCP.sh and compares the results between 
% Baseline and Baseline_GSR pipelines.
% 
% Inputs:
%   - LRR_dir
%     The output directory of CBIG_LiGSR_LRR_uniittest_PMAT_HCP.sh. It is
%     assumed that there are two subfolders call "Baseline" and "GSR" 
%     storing the kernel ridge regression results generated using each
%     preprocessing pipeline.
% 
% Outputs:
%     A file called [LRR_dir '/compare_2pipe/final_result.mat'] will be
%     created to store the integrated result and the comparison between the
%     two pipelines. It contains:
%     (1) mean_acc_dif: a 3x1 vector of fluid intelligence prediction
%                       accuracy difference between the two preprocessing
%                       pipelines, averaged within each random data split.
%     (2) acc_GSR_mean: a 3x1 vector of accuracy of fluid intelligence
%                       with global signal regression, averaged within each
%                       random data split.
%     (3) acc_Baseline_mean: a 3x1 vector of accuracy of fluid
%                            intelligence with baseline fMRI preprocessing,
%                            averaged within each random data split.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~exist(fullfile(LRR_dir, 'compare_2pipe'), 'dir'))
    mkdir(fullfile(LRR_dir, 'compare_2pipe'))
end

behaviors = {'PMAT24_A_CR'};
pipelines = {'GSR', 'Baseline'};

for b = 1:length(behaviors)
    for pipe = 1:length(pipelines)
        for seed = 1:3
            load(fullfile(LRR_dir, pipelines{pipe}, ['randseed_' num2str(seed)], ...
                'results', 'optimal_acc', [behaviors{b} '.mat']))
            
            if(pipe == 1)
                acc_GSR_mean(seed, b) = mean(acc_corr_test);
            else
                acc_Baseline_mean(seed, b) = mean(acc_corr_test);
            end
        end
    end
end
mean_acc_dif = acc_GSR_mean - acc_Baseline_mean;
save(fullfile(LRR_dir, 'compare_2pipe', 'final_result.mat'), 'mean_acc_dif', 'acc_GSR_mean', ...
    'acc_Baseline_mean')

end

