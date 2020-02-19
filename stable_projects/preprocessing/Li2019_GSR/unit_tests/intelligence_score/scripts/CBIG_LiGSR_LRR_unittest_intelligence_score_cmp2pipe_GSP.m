function CBIG_LiGSR_LRR_unittest_intelligence_score_cmp2pipe_GSP( LRR_dir )

% CBIG_LiGSR_LRR_unittest_intelligence_score_cmp2pipe_GSP( LRR_dir )
% 
% This is step 2 of the unit test which predicts two intelligences scores
% in the GSP dataset.
% 
% This function integrates the output files of
% CBIG_LiGSR_LRR_unittest_intelligence_score_cmp2pipe_GSP.sh and compares
% the results between Baseline and Baseline_GSR pipelines.
% 
% Inputs:
%   - KRR_dir
%     The output directory of CBIG_LiGSR_LRR_uniittest_intelligence_score_GSP.sh. 
%     It is assumed that there are two subfolders call "Baseline" and "GSR"
%     storing the linear ridge regression results generated using each
%     preprocessing pipeline.
% 
% Outputs:
%     A file called [LRR_dir '/compare_2pipe/final_result.mat'] will be
%     created to store the integrated result and the comparison between the
%     two pipelines. It contains:
%     (1) mean_acc_dif: a 3x2 matrix of the prediction accuracy difference 
%                       between the two preprocessing pipelines of Shipley
%                       Vocabulary and WAIS - Matrix Reasoning, averaged
%                       within each random data split.
%     (2) acc_GSR_mean: a 3x2 matrix of accuracy of Shipley Vocabulary and
%                       WAIS - Matrix Reasoning with global signal
%                       regression, averaged within each random data split.
%     (3) acc_Baseline_mean: a 3x2 matrix of accuracy of Shipley Vocabulary
%                            and WAIS - Matrix Reasoning with baseline fMRI
%                            preprocessing, averaged within each random
%                            data split.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

behaviors = {'Shipley_Vocab_Raw', 'Matrix_WAIS'};
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
mkdir(fullfile(LRR_dir, 'compare_2pipe'))
save(fullfile(LRR_dir, 'compare_2pipe', 'final_result.mat'), 'mean_acc_dif', 'acc_GSR_mean', ...
    'acc_Baseline_mean')

end

