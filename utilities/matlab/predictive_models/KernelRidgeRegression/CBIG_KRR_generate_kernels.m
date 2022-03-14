function CBIG_KRR_generate_kernels( feature_mat, sub_fold, outdir, ker_param, similarity_mat, cov_X )

% CBIG_KRR_generate_kernels( feature_mat, sub_fold, outdir, ker_param, similarity_mat, cov_X )
% 
% This function calculates the inner-loop and training-test kernels for
% kernel ridge regression. 
% 
% Inputs:
%   - feature_mat
%     Feature matrix.
%     Generally, "feature_mat" is a #features x #subjects matrix.
%     However if the user input a 3-D matrix, this function will assume
%     that the feature matrix is the connectivity matrix (eg. functional
%     connectivity across ROIs). Specifically, "feature_mat" will be a
%     #ROIs1 x #ROIs2 x #subjects matrix.
% 
%   - sub_fold
%     The data split for cross-validation.
%     It is a num_test_folds x 1 structure with a field "fold_index".
%     sub_fold(i).fold_index is a #subjects x 1 binary vector, where 1
%     refers to the test subjects in the i-th fold; 0 refers to the
%     training subjects in the i-th fold.
%     If the user does not need cross-validation, the length of sub_fold
%     (i.e. the number of folds) will be 1. In this case,
%     sub_fold.fold_index has 3 unique values:
%     0 - training set;
%     1 - test set;
%     2 - validation set.
% 
%   - outdir
%     A string, the full path of output directory.
%     For cross-validation case, 
%     (1) a subfolder [outdir '/FSM_innerloop'] will be created to store
%     inner-loop CV kernels.
%     (2) a subfolder [outdir '/FSM_test'] will be created to store
%     test CV kernels.
%     For training-validation-test case,
%     (1) a subfolder [outdir '/FSM_train'] will be created to store
%     training kernels.
%     (2) a subfolder [outdir '/FSM_validation'] will be created to store
%     validation kernels.
%     (3) a subfolder [outidr '/FSM_test'] will be created to store test
%     kernels.
% 
%   - ker_param (optional)
%     A K x 1 structure with two fields: type and scale. K denotes the
%     number of kernels.
%     ker_param(k).type is a string of the type of k-th kernel. Choose from
%                       'corr'        - Pearson's correlation;
%                       'Gaussian'    - Gaussian kernel;
%                       'Exponential' - exponential kernel.
%     ker_param(k).scale is a scalar specifying the scaling factor of k-th
%     kernel (only for Gaussian kernel or exponential kernel). 
%     If ker_param(k).type == 'corr', ker_param(k).scale = NaN.
%     
%     If "ker_param" is not passed in, only correlation kernel will be
%     calculated.
% 
%   - similarity_mat (optional)
%     A #subjects x #subjects inter-subject similarity matrix pre-computed
%     by the user. If this argument is passed in, this function will only
%     split this similarity matrix into kernels for different folds without
%     further operations. By passing in this argument, "feature_mat" is not
%     needed (you can pass in an empty matrix for the "feature_mat" input).
%
%   - cov_X (optional)
%     A #subjects x #covariates matrix, the covariates to be regressed from 
%     the features. The regression coefficients will be estimated from only
%     training set. The estimated coefficients then will be applied to the
%     test set.
% 
% Outputs:
%     Case 1. (cross-validation)
%     For each fold and each kernel hyperparameter, a #TrainingSubjects x
%     #TrainingSubjects matrix will be saved in [outdir '/FSM_innerloop']
%     folder as the inner-loop cross-validation kernel. 
%     Meanwhile, a #AllSubjects x #AllSubjects kernel matrix will be saved
%     in [outdir '/FSM_test'] folder for the training-test
%     cross-validation. (The ordering of subjects follows the original
%     subject list, i.e. the ordering in "feature_mat" or "similarity_mat".)
% 
%     Case 2 (training, validation, and test)
%     In this case, cross-validation is not performed. Instead, the data
%     are split into training, validation and test sets. For each kernel
%     hyperparameter, a (#TrainingSubjects + #ValidationSubjects) x
%     (#TrainingSubject + #ValidationSubjects) kernel matrix will be saved
%     in [outdir '/FSM_innerloop'] folder for the training and validation
%     phase. (Training subjects at first, then followed by validation
%     subjects.)
%     Meanwhile, a (#TrainingSubjects + #TestSubjects) x (#TrainingSubjects
%     + #TestSubjects) kernel matrix will be saved in [outdir '/FSM_test']
%     folder for the test phase. (Training subjects at first, then followed
%     by test subjects.)
% 
% Written by Jingwei, Ru(by) and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% If feature_mat is a connectivity matrix, extract the lower triangular
% entries and construct feature_mat as a #features x #subjects matrix.
if(size(feature_mat, 3) > 1)
    tril_ind = find(tril(ones(size(feature_mat,1), size(feature_mat,2)), -1) == 1);
    feature_mat = reshape(feature_mat, size(feature_mat,1)*size(feature_mat,2), size(feature_mat,3));
    feature_mat = feature_mat(tril_ind,:);
end

if(~exist('ker_param', 'var'))
    ker_param.type = 'corr';
    ker_param.scale = NaN;
end
num_ker = length(ker_param);

for k = 1:num_ker
    if(strcmp(ker_param(k).type, 'corr'))
        outstem{k} = '';
    else
        outstem{k} = ['_' num2str(ker_param(k).scale)];
    end
end

num_test_folds = length(sub_fold);
if(num_test_folds == 1)
    % training-validation-test case
    train_ind = sub_fold.fold_index == 0;
    test_ind = sub_fold.fold_index == 1;
    valid_ind = sub_fold.fold_index == 2;

    % regress covariates from features (if necessary)
    curr_outdir = fullfile(outdir, 'FSM_innerloop');
    if(~exist(curr_outdir, 'dir')); mkdir(curr_outdir); end
    if(exist('cov_X', 'var') && ~isempty(cov_X))
        fprintf('Regress covariates from features\n')
        cov_X_train_mean = mean(cov_X(train_ind, :));
        [feature_resid(train_ind, :), beta] = CBIG_regress_X_from_y_train(...
            feature_mat(:, train_ind)', cov_X(train_ind, :));
        feature_resid(valid_ind, :) = CBIG_regress_X_from_y_test(...
            feature_mat(:, valid_ind)', cov_X(valid_ind, :), beta, cov_X_train_mean);
        feature_resid(test_ind, :) = CBIG_regress_X_from_y_test(...
            feature_mat(:, test_ind)', cov_X(test_ind, :), beta, cov_X_train_mean);

        feature_mat = feature_resid';
        clear feature_resid
        save(fullfile(outdir, 'FSM_innerloop', 'beta.mat'), 'beta')
    end
    
    % training & validation
    fprintf('Training & validation FSM\n')
    for k = 1:num_ker
        fprintf('-- Current kernel type is %s, scale %f\n', ker_param(k).type, ker_param(k).scale);
        curr_outname = fullfile(curr_outdir, ['FSM_' ker_param(k).type outstem{k} '.mat']);
        if(~exist(curr_outname, 'file'))
            if(~exist('similarity_mat', 'var') || isempty(similarity_mat))
                FSM = CBIG_crossvalid_kernel_with_scale(feature_mat(:, train_ind), feature_mat(:, valid_ind), ...
                    [], [], ker_param(k).type, ker_param(k).scale);
            else
                FSM = similarity_mat(or(train_ind, valid_ind), or(train_ind, valid_ind));
                % use "+" because training and validation sets are supposed
                % to be exclusive.
            end
            save(curr_outname, 'FSM'); clear FSM
        else
            fprintf('Already exist. Skipping ...\n')
        end
    end
    
    % test
    fprintf('Test FSM\n')
    curr_outdir = fullfile(outdir, 'FSM_test');
    if(~exist(curr_outdir, 'dir')); mkdir(curr_outdir); end
    for k = 1:num_ker
        fprintf('-- Current kernel type is %s, scale %f\n', ker_param(k).type, ker_param(k).scale);
        curr_outname = fullfile(curr_outdir, ['FSM_' ker_param(k).type outstem{k} '.mat']);
        if(~exist(curr_outname, 'file'))
            if(~exist('similarity_mat', 'var') || isempty(similarity_mat))
                FSM = CBIG_crossvalid_kernel_with_scale(feature_mat(:, train_ind), feature_mat(:, test_ind), ...
                    [], [], ker_param(k).type, ker_param(k).scale);
            else
                FSM = similarity_mat(or(train_ind, test_ind), or(train_ind, test_ind));
                % use "+" because training and test sets are supposed 
                % to be exclusive.
            end
            save(curr_outname, 'FSM'); clear FSM
        else
            fprintf('Already exist. Skipping ...\n')
        end
    end
else
    % cross-validation case
    for test_fold = 1:num_test_folds
        train_ind = sub_fold(test_fold).fold_index == 0;
        test_ind = sub_fold(test_fold).fold_index == 1;
        
        % Regress covariates from features
        curr_outdir = fullfile(outdir, 'FSM_innerloop', ['fold_' num2str(test_fold)]);
        if(~exist(curr_outdir, 'dir')); mkdir(curr_outdir); end
        if(exist('cov_X', 'var') && ~isempty(cov_X))
            fprintf('Fold %d, regress covariates from features\n', test_fold)
            cov_X_train_mean = mean(cov_X(train_ind, :));
            [feature_resid(train_ind, :), beta] = CBIG_regress_X_from_y_train(...
                feature_mat(:, train_ind)', cov_X(train_ind, :));
            feature_resid(test_ind, :) = CBIG_regress_X_from_y_test(...
                feature_mat(:, test_ind)', cov_X(test_ind, :), beta, cov_X_train_mean);

            feature_mat = feature_resid';
            clear feature_resid
            save(fullfile(curr_outdir, 'beta.mat'), 'beta')
        end

        % inner-loop kernels
        fprintf('Fold %d, inner-loop FSM\n', test_fold);
        for k = 1:num_ker
            fprintf('-- Current kernel type is %s, scale %f\n', ker_param(k).type, ker_param(k).scale);
            curr_outname = fullfile(curr_outdir, ['FSM_' ker_param(k).type outstem{k} '.mat']);
            if(~exist(curr_outname, 'file'))
                if(~exist('similarity_mat', 'var') || isempty(similarity_mat))
                    FSM = CBIG_crossvalid_kernel_with_scale(feature_mat(:, train_ind), [], [], [], ...
                        ker_param(k).type, ker_param(k).scale);
                else
                    FSM = similarity_mat(train_ind, train_ind);
                end
                save(curr_outname, 'FSM'); clear FSM
            else
                fprintf('Already exist. Skipping ...\n')
            end
        end
        
        % traing-test kernels
        fprintf('Fold %d, test FSM\n', test_fold);
        curr_outdir = fullfile(outdir, 'FSM_test', ['fold_' num2str(test_fold)]);
        if(~exist(curr_outdir, 'dir')); mkdir(curr_outdir); end
        for k = 1:num_ker
            fprintf('-- Current kernel type is %s, scale %f\n', ker_param(k).type, ker_param(k).scale);
            curr_outname = fullfile(curr_outdir, ['FSM_' ker_param(k).type outstem{k} '.mat']);
            if(~exist(curr_outname, 'file'))
                if(~exist('similarity_mat', 'var') || isempty(similarity_mat))
                    FSM = CBIG_crossvalid_kernel_with_scale(feature_mat(:, train_ind), feature_mat(:, test_ind), ...
                        train_ind, test_ind, ker_param(k).type, ker_param(k).scale);
                else
                    train_test_ind = or(train_ind, test_ind);
                    FSM = similarity_mat(train_test_ind, train_test_ind);
                end
                save(curr_outname, 'FSM'); clear FSM
            else
                fprintf('Already exist. Skipping ...\n')
            end
        end
    end
end

end


