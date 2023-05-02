function CBIG_ICCW_weights_similarity(result_dir, output_dir)

% function CBIG_ICCW_haufe_similarity(result_dir, output_dir)
%
% This function computes the similarity between the original regression weights across
% regression models.
%
% Inputs:
%   - result_dir
%     This refers to the directory in which the results were saved. Assumes a sub-directory
%     for each regression model (KRR, LRR, LASSO and RF).
%
%   - output_dir
%     This refers to the directory to save the collated results.
%
% Outputs:
%   - sim_behav
%     A 3x3 matrix containing the similarity values between KRR, LRR, LASSO.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% setting up
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'ChenOoi2023_ICCW', 'analysis', 'utilities');
addpath(genpath(project_code_dir));

sim_behav = zeros(3,3,36);

% iterate over behaviors
for i = 1:36
    sim = zeros(3,3,252);
    % load KRR weights
    load(fullfile(result_dir, 'KRR', num2str(5260), 'weights' strcat('weights_score', num2str(i),'_all_folds.mat')));
    % iterate over folds
    for j = 1:252
        % extract weights for KRR, LRR and LASSO
        KRR_weights = weights_all_folds(:,j);
        load(fullfile(result_dir, 'LRR', num2str(5260), strcat(behav_' num2str(j)), 'rng1', 'params', ...
            strcat('fold_', num2str(i)), 'acc_test_all_score.mat'));
        LRR_weights = beta;
        load(fullfile(result_dir, 'LASSO', num2str(5260), strcat(behav_' num2str(j)), 'rng1', 'params', ...
            strcat('fold_', num2str(i)), 'acc_test_all_score.mat'));
        LASSO_weights = beta;
        
        % compare to KRR
        sim(2,1,j) = CBIG_ICCW_ICC_1_1([KRR_weights LRR_weights]);
        sim(3,1,j) = CBIG_ICCW_ICC_1_1([KRR_weights LASSO_weights]);
        % compare to LRR
        sim(3,2,j) = CBIG_ICCW_ICC_1_1([LRR_weights LASSO_weights]);
    end
    sim_behav(:,:,i) = mean(sim,3);
end
% save final results
save(fullfile(output_dir, 'sim_behav_weights.mat'), 'sim_behav')

rmpath(genpath(project_code_dir));
end