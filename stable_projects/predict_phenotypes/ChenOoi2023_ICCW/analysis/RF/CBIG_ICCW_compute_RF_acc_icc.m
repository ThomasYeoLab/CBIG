function CBIG_ICCW_compute_RF_acc_icc(result_dir, output_dir)

% function CBIG_ICCW_compute_RF_acc_icc(result_dir, output_dir)
%
% This function collates the accuracy results, and computes the ICC for the Haufe
% transformed weights, as well as the original regression weights for RF.
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
%   - acc_<sample_size>_RF
%     A mat file with a matrix of #folds/2 x #behaviours. Accuracy is in terms of
%     correlation averaged over the two split halves.
%
%   - icc_<sample_size>_RF_Haufe
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances (Haufe) values from the two split halves.
%
%   - icc_<sample_size>_RF_weights
%     A mat file with a matrix of #folds/2 x #behaviours. ICC is calculated
%     from the feature importances (regression weights) values from the two split halves.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% setting up
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'predict_phenotypes', 'ChenOoi2023_ICCW', 'analysis', 'utilities');
addpath(genpath(project_code_dir));

num_behav = 1; % change to 39 for all behaviors;
behavs = [1:num_behav];
s_sizes = [800 400]; % change to [5260 3000 2000 800 400] for all samples
RF_dir = fullfile(result_dir,'RF');
mkdirp(output_dir);

% save accuracies, Haufe transform and original weights ICC in mat files
for s = 1:length(s_sizes)
    disp(['Processing for sample size ', num2str(s_sizes(s))])
    N = num2str(s_sizes(s));
    % load fold information
    load(fullfile(result_dir, 'KRR', N, 'no_relative_5_fold_sub_list.mat'))
    icc_Haufe = zeros(126,num_behav);
    icc_weights = zeros(126,num_behav);
    acc = zeros(126,num_behav);
    
    % iterate over folds
    for i = 1:126
        fold = strcat('fold_', num2str(i));
        matching_fold = strcat('fold_', num2str(253-i));
        % iterate over behaviors
        for tmp = 1:num_behav
            j = behavs(tmp);
            behav = strcat('behav_', num2str(j));
            
            % get details for first split-half fold
            % get accuracies
            curr_path = fullfile(RF_dir, N, behav, 'rng0', fold);
            acc_csv = csvread(fullfile(curr_path, strcat('ABCD_RF_', behav, '_predacc.csv')),1,0);
            acc1 = acc_csv(1);
            % get Haufe transformed and condition variable importance values
            fi_csv = csvread(fullfile(curr_path, strcat('ABCD_RF_', behav, '_fi.csv')),1,0);
            beta1 = fi_csv(:,1);
            pfm1 = fi_csv(:,4);
            
            % get details for matching split-half fold
            % get accuracies
            matching_path = fullfile(RF_dir, N, behav, 'rng0', matching_fold);
            acc_csv = csvread(fullfile(matching_path, strcat('ABCD_RF_', behav, '_predacc.csv')),1,0);
            acc2 = acc_csv(1);
            % get Haufe transformed and condition variable importance values
            fi_csv = csvread(fullfile(matching_path, strcat('ABCD_RF_', behav, '_fi.csv')),1,0);
            beta2 = fi_csv(:,1);
            pfm2 = fi_csv(:,4);
            
            % get averaged accuracy and icc for regression weights and Haufe transform
            acc(i,j) = (acc1+acc2)/2;
            icc_weights(i,j) = CBIG_ICCW_ICC_1_1([beta1 beta2]);
            icc_Haufe(i,j) = CBIG_ICCW_ICC_1_1([pfm1 pfm2]);
        end
    end
    
    % save results
    save(fullfile(output_dir, strcat('acc_', N, '_RF.mat')), 'acc');
    icc = icc_Haufe;
    save(fullfile(output_dir, strcat('icc_', N, '_RF_Haufe.mat')), 'icc');
    icc = icc_weights;
    save(fullfile(output_dir, strcat('icc_', N, '_RF_weights.mat')), 'icc');
end

rmpath(genpath(project_code_dir));
end
