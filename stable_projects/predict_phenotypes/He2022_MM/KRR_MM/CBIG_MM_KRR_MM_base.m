function CBIG_MM_KRR_MM_base(CBIG_CODE_DIR, subject_list, phe_csv, FC_file, ...
    output_dir, rngs, phe_list, subj_meta_train, subj_meta_test)

% CBIG_MM_KRR_MM_base(CBIG_CODE_DIR, subject_list, phe_csv, FC_file, ...
%   output_dir, rngs, phe_list, subj_meta_train, subj_meta_test)
% 
% This function runs the kernel ridge regression algorithm to all of
% training meta-set non-brain-imaging phenotypes with training meta-set
% subjects.
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
% 
%   - subject_list
%     Full path of subject list of UK Biobank dataset. It should be a txt
%     file that contains #subject of line, while each line is the subject
%     id of 1 subject.
% 
%   - phe_csv
%     Full path of the phes csv for UK Biobank dataset. It should contains the
%     phes of each subjects.
%  
%   - FC_file
%     Full path of the functional connectivity matrix. A matrix "corr_mat" is
%     assumed to be saved in this file. "corr_mat" should be a 3-D matrix 
%     with dimension of #ROI x #ROI x #subjects. Since it is a connectivity
%     matrix, only the lower-triangular off-diagonal entries will be
%     considered as features because the connectivity matrix is symmetric. The
%     order of "corr_mat" should follow the subject order of "subject_list" 
% 
%   - output_dir
%     Full path of the output directory.
% 
%   - rngs
%     Number (integer) of random number generator repeats of kernel ridge
%     regression training, validation and testing split. It can be number or
%     string.
% 
%   - phe_list
%     Full path of non-brain-imaging phenotypes list to perform the
%     kernel ridge regression. Or it can be '31-0.0' which is sex that we 
%     seperately run because it is categorical.
%
%   - subj_meta_train
%     Full path of subject list of training meta-set. It should be a txt
%     file that contains #subject of line, while each line is the subject
%     id of 1 subject.
%
%   - subj_meta_test
%     Full path of subject list of test meta-set. It should be a txt
%     file that contains #subject of line, while each line is the subject
%     id of 1 subject.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    %% add kernel regression code to path
    addpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
        'predictive_models', 'KernelRidgeRegression')));

    %% generate setup file
    param = load(fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing',...
        'Li2019_GSR', 'examples', 'ref_output', 'KernelRidgeRegression',...
        'setup_file.mat'));
    
    % check FC_file to set pre_fix
    if endsWith(FC_file, 'exp_train_test_pfc.mat')
        pre_fix = 'exp';
    else
        pre_fix = 'ukbb';
    end

    % get matrix of phe
    csv_files = {phe_csv};
    subject_header = 'eid';

    if strcmp(phe_list, '31-0.0')
        y_names = {phe_list};
        need_znorm = false;
    else
        need_znorm = true;
        fileID = fopen(phe_list);
        temp = textscan(fileID, '%s');
        y_names = temp{1}';
        fclose(fileID)

        if any(strcmp(y_names, '31-0.0'))
            y_names = y_names(~strcmp(y_names, '31-0.0'));
            disp('you also need to run 31-0.0 (sex) individually')
        end

    end

    y_types = cell(1, size(y_names, 2));
    y_types(:) = {'continuous'};
    outname = fullfile(output_dir, ['phe_measures.mat']);
    delimiter = ',';

    if exist(outname, 'file')
        tmp = load(outname);
        y = tmp.y;
    else
        y = CBIG_read_y_from_csv(csv_files, subject_header, y_names, y_types, ...
            subject_list, outname, delimiter);
    end

    param.y = y;

    % only use 1 as covariates
    param.covariates = []; % 'none';

    % load FC matrix
    temp = load(FC_file);
    param.feature_mat = temp.corr_mat;

    % number of inner folds
    param.num_inner_folds = 1;

    % do not save out kernel
	param.save_kernel = 0;

    % range of lambda for kernel regression
    param.lambda_set = [0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 ...
                        0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];

    % metric for tune net
    param.metric = 'none';

    rng_num = str2num(rngs);

    data_dir = fullfile(output_dir, [pre_fix '_rng_num_' num2str(rng_num)]);

    if ~exist(data_dir, 'dir')
        mkdir(data_dir)
    end

    % check result file, if already have, pass this rng
    if isfile(fullfile(data_dir, 'final_result.mat'))
        disp(['rng ' num2str(rng_num) ' already have final_result.mat'])
        return
    end

    % split data
    CBIG_MM_KRR_MM_split_subject(data_dir, subject_list, rng_num, 0.2, subj_meta_train, subj_meta_test, pre_fix)
    temp = load(fullfile(data_dir, [pre_fix '_subject_split.mat']));
    param.sub_fold = temp.sub_fold;

    if need_znorm
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
    fsm_del = fullfile(data_dir, 'innerloop_cv');
    rmdir(fsm_del, 's')

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
        y_tra = y(ind == 0);
    else
        y_tra = y(ind == 0, :);
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
