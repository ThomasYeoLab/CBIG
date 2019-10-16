function kernel = CBIG_crossvalid_kernel_with_scale(data_KbyS1, data_KbyS2, index_S1, index_S2, type, scale)

% kernel = CBIG_crossvalid_kernel_with_scale(data_KbyS1, data_KbyS2, index_S1, index_S2, type, scale)
% 
% This function calculates the kernel matrix between subjects from two sets 
% (e.g. training set and test set) for predictive models (e.g. kernel ridge
% regression).
% 
% Inputs:
%   - data_KbyS1
%     A #features x (# subjects from set 1) matrix. It should correspond to
%     the features of the training subjects.
% 
%   - data_KbyS2
%     A #features x (# subjects from set 2) matrix.
%     It can be empty (i.e. []), the kernel function will only apply on
%     "data_KbyS1".
%     If it is not empty, it can be the features of the test subjects (to
%     compute the kernel for training-test cross-validation), or it can be
%     the features of the validation subjects (when cross-validation is not
%     performed and the the kernel will be computed for the validation phase).
% 
%   - index_S1
%     A binary vector with length of (# total subjects). An entry with 1
%     represents that this subject is in the first set of subjects, meaning
%     there is a column in data_KbyS1 that corresponds to this subject,
%     i.e. all_subjects_list(index_S1) = subjects_set1.
%     If it is empty (i.e. []), then the full subject list ordering will be
%     assumed to have the form [subjects_set1, subjects_set2], i.e. the
%     first S1 subjects are assigned to set 1, and the last S2 subjects are
%     assigned to set 2.
% 
%   - index_S2
%     A binary vector with length of (# total subjects). An entry with 1
%     represents that this subject is in the second set of subjects,
%     meaning there is a column in data_KbyS2 that corresponds to this
%     subject,
%     i.e. all_subjects_list(index_S2) = subjects_set2.
%     It it is empty (i.e. []), then the full subject list ordering will be
%     assumed to have the form [subjects_set1, subjects_set2], i.e. the
%     first S1 subjects are assigned to set 1, and the last S2 subjects are
%     assigned to set 2.
%  
%     Note: the 1s in index_S1 and index_S2 should be exclusive and
%     complementary, i.e. the same entry in index_S1 and index_S2 cannot be
%     both 1 (or both 0).
% 
%   - type
%     A string that determines the mode of tabulating the kernel. Choose from:
%     'corr'        - Pearson's correlation;
%     'Gaussian'    - Gaussian kernel;
%     'Exponential' - exponential kernel.
%     If Gaussian kernel or exponential kernel is chosen, data_KbyS1 and
%     data_KbyS2 will be normalized by the mean and standard deviation of
%     data_KbyS1 (i.e. the training set).
% 
%   - scale (optional)
%     A scalar, the scaling factor of Gaussian kernel or exponential kernel.
% 
% Outputs:
%   - kernel
%     A #all_subjects x #all_subjects kernel matrix, where #all_subjects =
%     #subjects_set1 + #subjects_set2.
%     If index_S1 and index_S2 are passed in, then the ordering of subjects
%     in "kernel" follows:
%     1. all_subjects_list(index_S1) = subjects_set1;
%     2. all_subjects_list(index_S2) = subjects_set2.
%     If index_S1 and index_S2 are not passed in, then the ordering of
%     subjects in "kernel" are the same as the ordering of data_KbyS1 and
%     data_KbyS2 (S1 first, then followed by S2).
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Ru(by) Kong and Jingwei Li

if(~strcmp(type,'corr'))
        
    K = size(data_KbyS1, 1); % feature number, same for data_KbyS1 and data_KbyS2
    SD1 = std(data_KbyS1,0,2); % std of the data_KbyS1
    mu1 = mean(data_KbyS1,2); % mean of the data_KbyS1


    % znormalize data_KbyS1 by SD1 and mu1
    data_KbyS1 = bsxfun(@minus, data_KbyS1, mu1);
    data_KbyS1 = bsxfun(@rdivide, data_KbyS1, SD1);
    if(~isempty(data_KbyS2))
        if(size(data_KbyS1, 1) ~= size(data_KbyS2, 1))
            error('Number of features for data_KbyS1 and data_KbyS2 must be the same!')
        end
        % znormalize data_KbyS2 also by SD1 and mu1
        data_KbyS2 = bsxfun(@minus, data_KbyS2, mu1);
        data_KbyS2 = bsxfun(@rdivide, data_KbyS2, SD1);
    end
end
% concatenate data_KbyS1 and data_KbyS2
if(isempty(index_S1) && isempty(index_S2))
    data_KbyS = [data_KbyS1, data_KbyS2];
else
    data_KbyS = zeros(size(data_KbyS1,1),length(index_S1));
    data_KbyS(:,index_S1) = data_KbyS1;
    data_KbyS(:,index_S2) = data_KbyS2;
end

%data_KbyS = zscore(data_KbyS,0,2);
if(strcmp(type, 'corr'))
    kernel = CBIG_nancorr(data_KbyS);
elseif(strcmp(type,'Exponential'))
    kernel = exp(-1*scale*squareform(pdist(data_KbyS'))/K);
elseif(strcmp(type,'Gaussian'))
    kernel = exp(-1*scale*squareform(pdist(data_KbyS').^2)/K);
else
    error('Unknown kernel type: %s', type);
end

end

