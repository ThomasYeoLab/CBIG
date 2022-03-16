function CBIG_MM_KRR_classical(CBIG_CODE_DIR, subject_list, phe_csv, FC_file, output_dir, rngs, phe, k, flag_pred_train)

% CBIG_MM_KRR_classical(CBIG_CODE_DIR, subject_list, phe_csv, FC_file, output_dir, rngs, phe, k, flag_pred_train)
% 
% This function runs the kernel ridge regression algorithm multiple times to one of
% phenotypes (non-brain-imaging phenotypes) for the input dataset.
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
% 
%   - subject_list
%     Full path of subject list of the dataset. It should be a txt
%     file that contains #subject of line, while each line is the subject
%     id of 1 subject.
% 
%   - phe_csv
%     Full path of the phenotypes csv for the dataset. It should contains the
%     phenotypes of each subjects.
%  
%   - FC_file
%     Full path of the functional connectivity matrix. A matrix "corr_mat" is
%     assumed to be saved in this file. "corr_mat" should be a 3-D matrix 
%     with dimension of #ROI x #ROI x #subjects. Since it is a connectivity
%     matrix, only the lower-triangular off-diagonal entries will be
%     considered as features because the connectivity matrix is symmetric. The
%     order of "corr_mat" should follow the subject order of "subject_list" 
% 
%   - rngs
%     Number (integer) of random number generator repeats of kernel ridge
%     regression training, validation and testing split. It can be number or
%     string.
% 
%   - phe
%     phenotypes (non-brain-imaging phenotypes) in dataset to perform the
%     kernel ridge regression.
% 
%   - output_dir
%     Full path of the output directory.
%
%   - k
%     Number (integer in string) of the K for K shot (participants) learning
%
%   - flag_pred_train
%     flag for haufe transform to predict train data
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    if endsWith(FC_file, 'HCP_diff_roi_pfc.mat')
        % check FC_file to set pre_fix
        pre_fix = 'HCP';
    elseif endsWith(FC_file, 'exp_test_pfc.mat')
        pre_fix = 'exp';
    else
        pre_fix = 'ukbb';
    end
    %% add kernel regression code to path
    addpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
        'predictive_models', 'KernelRidgeRegression')));

    if(isempty(flag_pred_train)) 
        flag_pred_train = false;
    elseif flag_pred_train == '1'
        flag_pred_train = true;
    else
        flag_pred_train = false;
    end

    %% generate setup file
    param = load(fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing',...
        'Li2019_GSR', 'examples', 'ref_output', 'KernelRidgeRegression',...
        'setup_file.mat'));

    % get output directory
    data_output_dir = fullfile(output_dir, ['output_phe_' phe]);
    if ~exist(data_output_dir, 'dir')
       mkdir(data_output_dir)
    end

    % get matrix of phenotype
    csv_files = {phe_csv};
    subject_header = 'eid';
    y_names = {phe};
    y_types = cell(1, size(y_names, 2));
    y_types(:) = {'continuous'};
    outname = fullfile(data_output_dir, ['phe_measures_k_' k '.mat']);
    delimiter = ',';
    y = CBIG_read_y_from_csv(csv_files, subject_header, y_names, y_types, ...
        subject_list, outname, delimiter);
    param.y = y;

    % only use 1 as covariates
    param.covariates = []; % 'none';

    % load FC matrix
    temp = load(FC_file);
    param.feature_mat = temp.corr_mat;

    % number of inner folds
    param.num_inner_folds = 5;

    % do not save out kernel
    param.save_kernel = 0;

    % range of lambda for kernel regression
    param.lambda_set = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 ...
        0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];

    % metric for tune net
    param.metric = 'none';

    if ischar(rngs)
        rngs = str2num(rngs);
    end
    rng_offset = 0;

    % main loop
    for rng_num = 1:rngs

        data_dir = fullfile(data_output_dir, [pre_fix '_' phe '_k_' k '_rng_num_' num2str(rng_num)]);
        if ~exist(data_dir, 'dir')
           mkdir(data_dir)
        end

        % check result file, if already have, pass this rng
        if isfile(fullfile(data_dir, 'final_result.mat'))
            disp([phe ' rng ' num2str(rng_num) ' already have final_result.mat'])
            continue
        end
        
        % split data
        CBIG_MM_KRR_classical_split_subject(data_dir, subject_list, rng_num + rng_offset, str2num(k), param.y, pre_fix)
        temp = load(fullfile(data_dir, [pre_fix '_subject_split.mat']));
        thr_y_znorm = 1.0 - 1.0/param.num_inner_folds;
        y_train = y(temp.sub_fold.fold_index == 0);
        while size(unique(y_train), 1) == 1 || max(hist(y_train,unique(y_train)))/size(y_train, 1) >= thr_y_znorm
            disp(['rng ' num2str(rng_num + rng_offset) ' for ' phe ' has ' num2str(thr_y_znorm) ' same value'])
            rng_offset = rng_offset + 1;
            CBIG_MM_KRR_classical_split_subject(data_dir, subject_list,...
                rng_num + rng_offset, str2num(k), param.y, pre_fix)
            temp = load(fullfile(data_dir, [pre_fix '_subject_split.mat']));
            y_train = y(temp.sub_fold.fold_index == 0);
        end
        param.sub_fold = temp.sub_fold;
        if ~strcmp(phe, '31-0.0')
            param.y = znorm(y, temp.sub_fold.fold_index);
        end

        % output directory
        param.outdir = data_dir;
        mkdir(param.outdir)

        %% run kernel regression
        CBIG_KRR_workflow_LITE(param, 0);

        % delete intermediate files
        fsm_del = fullfile(data_dir, 'test_cv');
        rmdir(fsm_del, 's')

        if flag_pred_train
            kernel = CBIG_KRR_generate_kernels_LITE(param.feature_mat,...
                param.outdir, param.save_kernel, param.ker_param);
            test_fold = 1;
            tmp_stem = '';
            tmp_fold = param.sub_fold;

            CBIG_MM_KRR_test_pred(test_fold, tmp_fold, ...
                param.outdir, tmp_stem, param.with_bias, kernel, ...
                param.ker_param, param.lambda_set, 'none');
            tmpdir = fullfile(param.outdir, 'test_cv', ['fold_' num2str(test_fold)]);
            tmpfile = fullfile(tmpdir, ['acc' tmp_stem '.mat']);
            tmptgt = fullfile(param.outdir, ['acc' tmp_stem '.mat']);
            copyfile(tmpfile, tmptgt);
            fsm_del = fullfile(data_dir, 'test_cv');
            rmdir(fsm_del, 's')

            acc_save = load(tmptgt);
            res_save = load(fullfile(param.outdir, ['final_result.mat']));
            ind = find(param.lambda_set == res_save.optimal_lambda);
            acc_save.optimal_y_t = acc_save.y_t{ind}{1};
            acc_save.optimal_y_p = acc_save.y_p{ind}{1};
            acc_save.optimal_acc = acc_save.acc{ind};
            acc_save.optimal_pred_stats = acc_save.pred_stats{ind};
            save(tmptgt,'-struct','acc_save')
        end

        fsm_del = fullfile(data_dir, 'innerloop_cv');
        rmdir(fsm_del, 's')

    end
    %% remove kernel regression code to path
    rmpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
        'predictive_models', 'KernelRidgeRegression')));

end


function ret = znorm(y, ind)
% ret = znorm(y, ind)
% 
% This function z normalize the y based on the entry with ind index
%
% Inputs:
%   - y
%     array of y, the data to be z-normalized, it can be vector or 2D matrix
%
%   - ind
%     array of index for y to calculate mean and std for z normalization. 
%
% Outputs:
%	- ret = z normalized y
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
