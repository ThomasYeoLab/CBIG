function [y_resid_perfold] = CBIG_crossvalid_regress_covariates_from_y( ...
    y_in, regressors, sub_fold, outdir, outstem)

% [y_resid_perfold] = CBIG_crossvalid_regress_covariates_from_y( ...
%     y_in, regressors, sub_fold, outdir, outstem)
% 
% This function regress the "regressors" from the target variables "y_in"
% for cross-validation based model usage. The regression coefficients will
% be generated from each training fold specified by "sub_fold" structure.
% The regression coefficients will be applied to both training and test
% folds.
% 
% Inputs:
%   - y_in
%     A #subjects x #TargetVariables matrix, where the regressors will be
%     regressed from. Examples of y_in could be the behavioral measures of
%     all subjects.
% 
%   - regressors
%     A #subjects x #covariates matrix of regressors.
%     If the user does not want to regress anything except performing
%     demean, he/she can pass in empty matrix as "regressors". The scripts
%     will automatically regress a vector of ones from y.
%     If the user doea not want to regress anything including the mean,
%     he/she should pass in 'NONE' so that this function can directly
%     assign the residual to be the same as the input.
% 
%   - sub_fold
%     The data split for cross-validation.
%     It is a num_test_fold x 1 structure with a field "fold_index".
%     sub_fold(i).fold_index is a #subjects x 1 binary vector, where 1
%     refers to the test subjects in the i-th fold; 0 refers to the
%     training subjects in the i-th fold.
%     If num_test_fold = 1, that means cross-validation will not be
%     performed. The datasets are assumed to be divided into training,
%     validation, and test sets. In this case, the regression coefficients
%     will also be computed based on the training set and applied to the
%     test set.
% 
%   - outdir
%     A string, the full path of output directory. A subfolder 
%     [outdir '/y/fold_' test_fold] will be created for each test fold.
% 
%   - outstem (optional)
%     A string appended to filenames to indentify the output file. For
%     example, if the output filename is
%     <path_to_file>/y_regress_58behaviors.mat, the outstem =
%     '58behaviors'. If outstem is empty, the output filename will be
%     <path_to_file>/y_regress.mat
% 
% Outputs:
%   - y_resid_perfold
%     Cell arrays with length of num_test_fold. Each cell array is a
%     #subjects x #TargetVariables matrix of the target variables after
%     regression. Each cell array will also be saved as a separate .mat
%     file under the folder [outdir, '/y/fold_' num2str(test_fold)] (see
%     the description of outstem).
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Authors: Jingwei Li, Ru(by) Kong

y_orig = y_in; 
num_test_folds = length(sub_fold);
for test_fold = 1:num_test_folds
    curr_outdir = fullfile(outdir, 'y', ['fold_' num2str(test_fold)]);
    mkdir(curr_outdir)
    if( exist('outstem', 'var') && ~isempty(outstem))
        outname = fullfile(curr_outdir, ['y_regress_' outstem '.mat']);
    else
        outname = fullfile(curr_outdir, ['y_regress.mat']);
    end
    
    if(~exist(outname, 'file'))
        if(ischar(regressors) && strcmpi(regressors, 'none'))
            y_resid = y_orig;
            beta = [];
        else
            beta = zeros(size(regressors,2)+1, size(y_in,2));
            for i = 1:size(y_in, 2)
                train_ind = sub_fold(test_fold).fold_index==0;
                test_ind = sub_fold(test_fold).fold_index==1;
                
                if(isempty(regressors))
                    X_train = [];
                    X_test = [];
                    X_train_mean = [];
                else
                    X_train = regressors(train_ind,:);
                    X_train_mean = mean(X_train);
                    X_test = regressors(test_ind,:);
                end
                
                [y_resid(train_ind,i), beta(:,i)] = ...
                    CBIG_regress_X_from_y_train(y_in(train_ind,i), X_train);
                y_resid(test_ind,i) = ...
                    CBIG_regress_X_from_y_test(y_in(test_ind,i), X_test, beta(:,i), X_train_mean);
                
                if(num_test_folds==1)
                    valid_ind = sub_fold(test_fold).fold_index==2;
                    if(isempty(regressors))
                        X_valid = [];
                    else
                        X_valid = regressors(valid_ind,:);
                    end
                    y_resid(valid_ind,i) = CBIG_regress_X_from_y_test(y_in(valid_ind,i), X_valid, ...
                        beta(:,i), X_train_mean);
                end
            end
        end
        
        save(outname, 'y_resid', 'y_orig', 'beta');
    else
        fprintf('Already exist. Skipping ...\n')
        load(outname)
    end
    
    y_resid_perfold{test_fold} = y_resid;
    clear y_resid
end

end

