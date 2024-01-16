function CBIG_MMP_HCP_compute_stacking_perm_stats(single_dir, stacking_dir, outstem, feature_mats_cell, ...
    score_ind, perm_seed_start, N_perm, family_list, subject_list, outdir)

% function CBIG_MMP_HCP_compute_stacking_perm_stats(single_dir, stacking_dir, outstem, ...
%    feature_mats_cell, score_ind, perm_seed_start, N_perm, site_list, outdir)
%
% This function performs the permutation test for the stacking model and
% saves out the prediction accuracies of the permutation test
%
% Inputs:
%   - single_dir:
%     The directory wherre the singleKRR results are.
%
%   - stacking_dir:
%     The directory wherre the stacking results are
%
%   - outstem
%     A string appended to the filename to specify the output files in the
%     singleKRR prediction. Must be the same as the input argument `outstem`
%     when you perform singleKRR
%
%   - feature_mats_cell
%     Cell of base file name for the input features.
%
%   - score_ind
%     A scalar. The index of the score you want to perfom the permutation test.
%     Range from 1 to # Target Variables in your kernel regression
%
%   - perm_seed_start
%     A scalar. The function will do the permutation N_perm times using
%     permutation seed from perm_seed_start to perm_seed_start+N_perm-1
%
%   - N_perm
%     Number of permutation
%
%   - family_list
%     A csv containing the family IDs of subjects.
%
%   - subject_list
%     A text file where each line is the subject ID to be extracted from family_list.
%
%   - outdir
%     Output directory for the permutation results
%
% Outputs:
%    prediction accuracies of the permutation test will be saved in the
%    outdir
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% check input arguments
if(isstr(score_ind))
    score_ind = str2num(score_ind);
end

if ~exist(outdir,'dir')
    mkdir(outdir);
end

% get subjects
fid = fopen(subject_list,'r');
file = textscan(fid, '%d');
fclose(fid);
subjects = file{:};

% get families
sub_table = readtable(family_list);
[~, loc] = ismember(subjects,sub_table.Subject);
family_all = sub_table.Family_ID(loc);
[~, ~, family_ind] = unique(family_all);
N_family = max(family_ind);
for i = 1:N_family
    family_count(i,:) = sum(family_ind == i);
end

%% load cross-validation split
N_seeds = 60;
variable_name = strcat('variable_', num2str(score_ind));
lvl1_lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];
for h=1:N_seeds
    seed = strcat('seed_',num2str(h));
    disp(seed);
    if ~exist(fullfile(outdir,seed),'dir')
        mkdir(fullfile(outdir,seed));
    end
    % get first level fold for KRR rsfmri
    sub_fold_file = dir(fullfile(stacking_dir, '..', 'KRR_features_rs', seed, 'results', '*sub_list.mat'));
    if length(sub_fold_file) ~= 1
        error('There should be one and only one sub_fold file in the singleKRR directory');
    end
    load(fullfile(stacking_dir, '..', 'KRR_features_rs', seed, 'results', sub_fold_file.name));
    N_fold = length(sub_fold);
    
    %% compute permutation
    metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
    for i = 1:length(metrics)
        stats_perm.(metrics{i}) = zeros(N_fold,N_perm);
    end
    
    outfile = fullfile(outdir, seed,...
        ['/acc_score' num2str(score_ind) '_allFolds_permStart' num2str(perm_seed_start) '.mat']);
    if ~exist(outfile)
        for i = 1:N_fold
            disp(['fold ' num2str(i)]);
            fold_name = strcat('fold_', num2str(i));
            % load data
            load(fullfile(stacking_dir, variable_name, seed, fold_name, ...
                'params', 'fold_1', ['selected_parameters_' outstem '.mat']));
            load(fullfile(stacking_dir, variable_name, seed, fold_name, ...
                'y', 'fold_1', ['y_regress_' outstem '.mat']));
            %y_resid = y_resid(:,score_ind);
            opt_lambda = curr_lambda;
            test_ind = sub_fold(i).fold_index;
            train_ind = ~test_ind;
            
            % load FSM
            all_feat_preds = [];
            for n = 1:size(feature_mats_cell,2)
                feature_name = feature_mats_cell{n};
                fprintf('%s \n', feature_name);
                feature_outstem = strcat('KRR_',feature_name);
                feature_result_mat = strcat('acc_KRR_',feature_name,'.mat');
                feat_train_tmp = load(fullfile(single_dir,feature_outstem,seed, ...
                    'results','innerloop_cv',fold_name,feature_result_mat));
                feat_test = load(fullfile(single_dir,feature_outstem,seed, ...
                    'results','test_cv',fold_name,feature_result_mat));
                feat_opt = load(fullfile(single_dir,feature_outstem,seed, ...
                    'results', strcat('final_result_',feature_outstem,'.mat')));
                feat_lambda = feat_opt.optimal_lambda(i, score_ind);
                % unscramble feat_train - needs to be updated if KRR cv rng changes
                num_inner_folds = 10;
                rng(1, 'twister');
                cv_idx = cvpartition(sum(train_ind), 'kfold', num_inner_folds);
                scrambled_idx = 1;
                clear feat_train
                for fol = 1:num_inner_folds
                    last_idx = scrambled_idx + cv_idx.TestSize(fol) - 1;
                    feat_train(cv_idx.test(fol)) = ...
                        feat_train_tmp.y_pred{:, ...
                            find(lvl1_lambda_set == feat_lambda)}{:,score_ind}(scrambled_idx:last_idx);
                    scrambled_idx = last_idx + 1;
                end
                feat_pred_y(train_ind) = feat_train;
                feat_pred_y(test_ind) = feat_test.y_p{:,find(lvl1_lambda_set == feat_lambda)}{:,score_ind};
                all_feat_preds = [all_feat_preds; feat_pred_y];
            end
            all_feat_preds = all_feat_preds';
            
            feat_train = all_feat_preds(train_ind,:);
            feat_pred = all_feat_preds(test_ind,:);
            
            % number of subjects
            num_train_subj = size(feat_train,1);
            num_test_subj = size(feat_pred,1);
            
            % Perform training and prediction
            for j = 1:N_perm
                % permute behavior within site
                y_perm = y_resid;
                % split into max(family_count) groups
                perm_group = zeros(length(family_ind),1);
                rng(j+perm_seed_start-1);
                for k = 1:N_family
                    sorting_order = datasample([1:max(family_count)],family_count(k),'Replace', false);
                    perm_group(find(family_ind == k)) = sorting_order;
                end
                % permute within groups
                for n = 1:max(family_count)
                    y_tmp = y_perm(perm_group == n);
                    y_tmp = y_tmp(randperm(length(y_tmp)));
                    y_perm(perm_group == n) = y_tmp;
                end
                
                y_train_resid = y_perm(train_ind);
                y_test_resid = y_perm(test_ind);
                
                % train and predict
                [acc, optimal_stats,pred,pred_train] = ...
                    CBIG_LRR_fitrlinear_train_test( feat_train, feat_pred, y_train_resid, y_test_resid, ...
                    opt_lambda, 'predictive_COD' );
                for k = 1:length(metrics)
                    stats_perm.(metrics{k})(i,j) = optimal_stats.(metrics{k});
                end
            end
        end
        save(outfile,'stats_perm');
    end
end
