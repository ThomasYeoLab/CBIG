function kernel = CBIG_KRR_generate_kernels_LITE(feature_mat, outdir, save_kernel, ker_param, similarity_mat)

% CBIG_KRR_generate_kernels_LITE( feature_mat, sub_fold, outdir, ker_param, similarity_mat )
% 
% This function calculates the inner-loop and training-test kernels for
% kernel ridge regression. 
% 
% This is a lite version of CBIG_KRR_generate_kernels.m. In this version,
% the kernels are not computed and saved separately for each inner-loop
% fold and each outer training-test fold. Instead, the kernel is only
% computed once across all the subjects, and for each fold, the scripts
% will grab the corresponding entries on the fly.
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
%     A subfolder [outdir '/FSM'] will be created to store the kernels
%     across all subjects.
%
%   - save_kernel
%     A string (choose from '0' or '1') or a scalar (choose from 0 or 1).
%     - save_kernel = 0 (or '0') means the algorithm is do not save kernel
%     into files.
%     - save_kernel = 1 (or '1') means the algorithm is save kernel into
%     files.
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
% Outputs:
%   - kernel
%     a #AllSubjects x #AllSubjects kernel matrix (The ordering of subjects 
%     follows the original subject list, i.e. the ordering in "feature_mat"
%     or "similarity_mat".)
%     if input 'save_kernel' = 1 (or '1'), this kernel matrix will be also
%     saved in [outdir '/FSM'].
%     if it has already been saved in [outdir '/FSM'], then a 'none' will 
%     be passed out. Later algorithm will load it from the file.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

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
    fprintf('-- Current kernel type is %s, scale %f\n', ker_param(k).type, ker_param(k).scale);
    if(strcmp(ker_param(k).type, 'corr'))
        outstem{k} = '';
    else
        outstem{k} = ['_' num2str(ker_param(k).scale)];
    end
    outname = fullfile(outdir, 'FSM', ['FSM_' ker_param(k).type outstem{k} '.mat']);
    if(save_kernel~=0)
        if(~exist(fullfile(outdir, 'FSM')))
            mkdir(fullfile(outdir, 'FSM'));
        end
    end
    
    %% Depending on kernel parameters, compute kernels
    if(~exist(outname, 'file'))
        if(~exist('similarity_mat', 'var') || isempty(similarity_mat))
        
            if(strcmp(ker_param(k).type, 'corr'))
                FSM = CBIG_self_corr(feature_mat);
                if(~isempty(isnan(feature_mat)))
                    nan_set = find(isnan(sum(feature_mat,1)));
                    for nan_s = nan_set
                        for all_s = 1:size(feature_mat,2)
                            nan_mask = ~isnan(feature_mat(:,nan_s)) & ~isnan(feature_mat(:,all_s));
                            FSM(nan_s,all_s) = CBIG_corr(feature_mat(nan_mask,nan_s),feature_mat(nan_mask,all_s));
                            FSM(all_s,nan_s) = FSM(nan_s,all_s);
                        end
                    end
                end
            else
                error('Only correlation kernel is supported!');
            end
        else
            FSM = similarity_mat;
        end
        tmp.FSM = FSM;
        kernel{k} = tmp;
        if(save_kernel~=0)
            save(outname, 'FSM'); clear FSM
        end
    else
        kernel = 'none';
        fprintf('Already exist. Skipping ...\n')
    end
end

end

