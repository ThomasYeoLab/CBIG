function CBIG_MMP__HCP_compute_singleKRR_perm_stats(singleKRR_dir, outstem, score_ind, ...
    perm_seed_start, N_perm, family_list, subject_list, outdir)

% function CBIG_MMP_HCP_compute_singleKRR_perm_stats(singleKRR_dir, outstem, score_ind, ...
% perm_seed_start, N_perm, family_list, subject_list, outdir)
%
% Adapted from CBIG_TRBPC_compute_singleKRR_perm_stats.
% This function performs the permutation test for the singleKRR model and
% saves out the prediction accuracies of the permutation test
%
% Inputs:
%   - singleKRR_dir:
%     The directory wherre the singleKRR results are
%
%   - outstem
%     A string appended to the filename to specify the output files in the
%     singleKRR prediction. Must be the same as the input argument `outstem`
%     when you perform singleKRR
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

% check if output already exist
%result_file = fullfile(outdir,['/acc_score' num2str(score_ind) '_allFolds_permStart' num2str(perm_seed_start) '.mat']);
%if ~exist(result_file,'file')
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
for h=1:N_seeds
    seed = strcat('seed_',num2str(h));
    disp(seed);
    if ~exist(fullfile(outdir,seed),'dir')
        mkdir(fullfile(outdir,seed));
    end
    sub_fold_file = dir(fullfile(singleKRR_dir, seed, 'results', '*sub_list.mat'));
    if length(sub_fold_file) ~= 1
        error('There should be one and only one sub_fold file in the singleKRR directory');
    end
    load(fullfile(singleKRR_dir, seed, 'results', sub_fold_file.name));
    N_fold = length(sub_fold);
    
    %% compute permutation
    metrics = {'corr','COD','predictive_COD','MAE','MAE_norm','MSE','MSE_norm'};
    for i = 1:length(metrics)
        stats_perm.(metrics{i}) = zeros(N_fold,N_perm);
    end
    
    load(fullfile(singleKRR_dir, seed, 'results', ['final_result_' outstem '.mat']));
    outfile = fullfile(outdir, seed,...
        ['/acc_score' num2str(score_ind) '_allFolds_permStart' num2str(perm_seed_start) '.mat']);
    if ~exist(outfile)
        for i = 1:N_fold
            disp(['fold ' num2str(i)]);
            % load data
            load(fullfile(singleKRR_dir, seed, 'results', 'y', ['fold_' num2str(i)], ['y_regress_' outstem '.mat']));
            y_resid = y_resid(:,score_ind);
            opt_lambda = optimal_lambda(i,score_ind);
            test_ind = sub_fold(i).fold_index;
            train_ind = ~test_ind;
            
            % load FSM
            load(fullfile(singleKRR_dir,seed, 'results', 'FSM_test', ['fold_' num2str(i)], 'FSM_corr.mat'));
            
            FSM_train = FSM(train_ind,train_ind,:);
            FSM_pred = FSM(test_ind,train_ind,:);
            
            % number of subjects
            num_train_subj = size(FSM_train,1);
            num_test_subj = size(FSM_pred,1);
            
            % Perform training and prediction
            K_r = FSM_train + opt_lambda*eye(num_train_subj);
            % training
            X = ones(num_train_subj,1);
            inv_K_r = inv(K_r);
            beta_pre = (X' * (inv_K_r * X)) \ X' * inv_K_r;
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
                % predict
                beta = beta_pre * y_train_resid;
                alpha = inv_K_r * (y_train_resid - X * beta);
                
                y_predicted = FSM_pred * alpha + ones(num_test_subj,1) .* beta;
                for k = 1:length(metrics)
                    stats_perm.(metrics{k})(i,j) = ...
                        CBIG_compute_prediction_acc_and_loss(y_predicted,y_test_resid,metrics{k},y_train_resid);
                end
            end
        end
        save(outfile,'stats_perm');
    end
end
