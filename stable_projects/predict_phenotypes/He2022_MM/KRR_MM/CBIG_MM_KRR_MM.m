function [opt_entry, opt_k, opt_test_corr, opt_test_cod] = ...
    CBIG_MM_KRR_MM(data_dir, krr_classical_dir, y_pred_all, test_csv, phe_tra,...
        phe_tes, phe_tes_mat, rng_num, k, prefix, eval_func)

% [opt_entry, opt_k, opt_test_corr, opt_test_cod] = ...
%   CBIG_MM_KRR_MM(data_dir, y_pred_all, test_csv, phe_tra,...
%       phe_tes, phe_tes_mat, rng_num, k, prefix, eval_func)
% 
% This function runs the meta-matching with kernel ridge regression
% algorithm to all testing non-brain-imaging phenotypes.
%
% Inputs:
%   - data_dir
%     Path of the your output and intermediate values. You can
%     also change this to any place you want.
%
%   - krr_classical_dir
%     Path of the the KRR classcial output and intermediate values.
%
%   - y_pred_all
%     Matrix with shape of #test meta-set subject * #training meta-set phenotypes.
%     It contains the prediction on training meta-set phenotypes of test meta-set
%     subjects by the kernel ridge regression model trained on training meta-set.
% 
%   - test_csv
%     Matlab table of the phenotpyes csv for test meta-set. It should contains
%     the phenotpyes of each subjects.
%  
%   - phe_tra
%     Cell (string) of training meta-set phenotypes names.
% 
%   - phe_tes
%     Cell (string) of test meta-set phenotypes names.
% 
%   - phe_tes_mat
%     Cell (string) of test meta-set phenotypes names for matlab table 
%     extraction.
% 
%   - rng_num
%     Number (integer) of random number generator to load from classical 
%     kernel ridge regression experiment.
%
%   - k
%     Number (integer) of the K for K shot (participants) learning
%
%   - prefix
%     str of the prefix for the dataset.
%
%   - eval_func
%     String for evaluation method for Meta-matching, can be 'corr' for 
%     correlation and 'cod' for COD.
%
% Outputs:
%   - opt_entry
%     A cell array with dimension of #test meta-set phenotypes. It contains
%     the optimal (matched) training meta-set phenotype for each test
%     meta-set phenotypes.
%
%   - opt_k
%     A cell array with dimension of #test meta-set phenotypes. It contains
%     the optimal (matched) training meta-set phenotype prediction for each
%     test meta-set phenotypes on k subjects.
%
%   - opt_test_corr
%     A cell array with dimension of #test meta-set phenotypes. It contains
%     the optimal (matched) test meta-set phenotype prediction correlation 
%     for each test meta-set phenotypes on remaining test subjects from test
%     meta-set (exclude k subjects).
%
%   - opt_test_cod
%     A cell array with dimension of #test meta-set phenotypes. It contains
%     the optimal (matched) test meta-set phenotype prediction COD for each
%     test meta-set phenotypes on remaining test subjects from test meta-set
%     (exclude k subjects).
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if strcmp(eval_func, 'corr')
        eval_fun = @corr;
    elseif strcmp(eval_func, 'cod')
        eval_fun = @cod;
    else
        error('unseen eval_func');
    end

    rng(rng_num, 'twister');

    opt_test_corr = zeros(length(phe_tes), 1);
    opt_test_cod = zeros(length(phe_tes), 1);
    opt_k = ones(length(phe_tes), 1) .* -Inf;
    for i = 1:length(phe_tes)

        y_test = test_csv{:, phe_tes_mat{i}};
        split_file_folder = fullfile(krr_classical_dir, ['output_phe_', phe_tes{i}], ...
            [prefix '_', phe_tes{i}, '_k_', num2str(k), '_rng_num_', num2str(rng_num)]);
        split_file = fullfile(split_file_folder, [prefix '_subject_split.mat']);
        k_split = load(split_file);
        k_index = logical((k_split.sub_fold.fold_index == 0)  + (k_split.sub_fold.fold_index == 2));
        test_index = (k_split.sub_fold.fold_index == 1);

        % z norm y_test
        y_test = znorm(y_test, ~k_index);

        for j = 1:length(phe_tra)
            y_pred = y_pred_all(:, j);

            tmp = eval_fun(y_pred(k_index), y_test(k_index));
            y_pred_tmp = -y_pred;
            tmp1 = eval_fun(y_pred_tmp(k_index), y_test(k_index));
            if tmp1 > tmp
                tmp = tmp1;
                y_pred = y_pred_tmp;
            end
            if tmp >= opt_k(i)
                opt_entry{i} = phe_tra{j};
                opt_k(i) = tmp;
                y_pred_opt = y_pred;
            end
        end
        if ~exist('opt_entry','var')
            disp('test')
        end
        disp([phe_tes{i}, ' ', opt_entry{i}, ' ', num2str(opt_k(i))])
        y_test_remain = y_test(test_index);
        real_index = ~isnan(y_test_remain);
        y_test_remain = y_test_remain(real_index);
        y_pred_remain = y_pred_opt(test_index);
        y_pred_remain = y_pred_remain(real_index);
        opt_test_corr(i) = corr(y_pred_remain, y_test_remain);
        opt_test_cod(i) = cod(y_pred_remain, y_test_remain);
        disp([num2str(opt_test_cod(i)), ' ', num2str(opt_test_corr(i)), ' ', num2str(tmp)])
    end

    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    res_dir = fullfile(data_dir, 'meta_result');
    mkdir(res_dir)
    save(fullfile(res_dir, ['final_result_with_' num2str(k) '_shot.mat']),...
        'phe_tes', 'opt_entry', 'opt_k', 'opt_test_corr', 'opt_test_cod')
end


function acc_out = reorganize_acc(acc)
    for i = 1:size(acc,1)
        for j = 1:size(acc,2)
            acc_each = acc{i,j};
            for m = 1:length(acc_each)
                acc_out(i,j,m) = acc_each(m);
            end
        end
    end
end


function r_square = cod(x, y)
    ss_res = sum((x - y).^2);
    ss_tot = sum((y).^2);
    r_square = 1 - ss_res/ss_tot;
end


function ret = znorm(y, ind)
    if isvector(y)
        y_tra = y(ind==0);
    else
        y_tra = y(ind==0, :);
    end
    y_mu = nanmean(y_tra, 1);
    y_std = nanstd(y_tra, 1);
    
    ret = bsxfun(@minus, y, y_mu);
    if sum(~y_std)
        warning('zero in standard deviation')
    else
        ret = bsxfun(@rdivide, ret, y_std);
    end
end