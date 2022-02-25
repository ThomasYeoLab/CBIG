function CBIG_MM_KRR_filter(CBIG_CODE_DIR, subject_list, phe_csv, FC_file, rngs, phe, output_dir)

% CBIG_MM_KRR_filter(CBIG_CODE_DIR, subject_list, phe_csv, FC_file, output_dir, rngs, phe)
% 
% This function runs the kernel ridge regression algorithm multiple times to one of
% phenotypes (non-brain-imaging phenotypes) for UK Biobank dataset for filter
% phenotypes. 
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
%     Full path of the phenotypes csv for UK Biobank dataset. It should contains the
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
%     phe (non-brain-imaging phenotypes) in UK Biobank dataset to perform the
%     kernel ridge regression.
% 
%   - output_dir
%     Full path of the output directory.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    %% add kernel regression code to path
    addpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
        'predictive_models', 'KernelRidgeRegression')));

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
    outname = fullfile(data_output_dir, 'phe_measures.mat');
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
    param.num_inner_folds = 1;

    % do not save out kernel
    param.save_kernel = 0;

    % range of lambda for kernel regression
    param.lambda_set = [0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 ...
        0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];

    % metric for tune net
    param.metric = 'none';

    if ischar(rngs)
        rngs = str2num(rngs);
    end

    % main loop
    for rng_num = 1:rngs

        % get the output directory for this rng
        data_dir = fullfile(data_output_dir, ['ukbb_' phe '_rng_num_' num2str(rng_num)]);
        if ~exist(data_dir, 'dir')
           mkdir(data_dir)
        end

        % check result file, if already have, pass this rng
        if isfile(fullfile(data_dir, 'final_result.mat'))
            disp([phe ' rng ' num2str(rng_num) ' already have final_result.mat'])
            continue
        end
        
        % split data
        CBIG_MM_KRR_filter_split_subject(data_dir, subject_list, rng_num, 0.2)
        temp = load(fullfile(data_dir, 'ukbb_subject_split.mat'));
        param.sub_fold = temp.sub_fold;

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

    end
    %% remove kernel regression code to path
    rmpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
        'predictive_models', 'KernelRidgeRegression')));

end

function subj_list = read_sub_list(subject_text_list)
% this function will output a 1xN cell where N is the number of
% subjects in the text_list, each subject will be represented by one
% line in the text file
% NOTE: multiple runs of the same subject will still stay on the same
% line
% Each cell will contain the location of the subject, for e.g.
% '<full_path>/subject1_run1_bold.nii.gz <full_path>/subject1_run2_bold.nii.gz'
    fid = fopen(subject_text_list, 'r');
    i = 0;
    while(1);
        tmp = fgetl(fid);
        if(tmp == -1)
            break
        else
            i = i + 1;
            subj_list{i} = tmp;
        end
    end
    fclose(fid);
end