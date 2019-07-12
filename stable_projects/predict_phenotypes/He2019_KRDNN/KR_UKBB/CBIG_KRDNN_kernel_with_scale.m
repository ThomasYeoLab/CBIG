function FSM=CBIG_KRDNN_kernel_with_scale(data_KbyS1,data_KbyS2,index_S1,index_S2,type,scale)
% CBIG_KRDNN_kernel_with_scale(data_KbyS1,data_KbyS2,index_S1,index_S2,type,scale)
% 
% This function compute kernel for kernel ridge regression algorithm based on
% data from two sets.
%
% Inputs:
%   - data_KbyS1
%     data of set 1, it is a KxS1 matrix, where K is the number of features,
%     S1 is the number of subjects in set 1.
%   - data_KbyS2
%     data of set 2, it is a KxS2 matrix, where K is the number of features,
%     S2 is the number of subjects in set 2. data_KbyS2 can be empty, then
%     Kernel function will only apply on data_KbyS1
%   - index_S1
%     index vector of set 1, it has S1 elements
%   - index_S2
%     index vector of set 2, it has S2 elements
%   - type
%     type of kernel to be generated. it has three options: 'corr'
%     (Pearson's correlation), 'Gaussian' (Gaussian kernel) and 'Exponential'
%     (exponential kernel). Default is 'corr'.
%   - scale
%     is a scalar specifying the scale (for Gaussian kernel or exponential
%     kernel). If type == 'corr', scale is not used.
%
% Outputs:
%   - FSM
%     The kernel calculated for kernel ridge regression.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~strcmp(type,'corr'))
        
    K = size(data_KbyS1, 1); % feature number, same for data_KbyS1 and data_KbyS2
    SD1 = std(data_KbyS1,0,2); % std of the data_KbyS1
    mu1 = mean(data_KbyS1,2); % mean of the data_KbyS1

    % znormalize data_KbyS1 by SD1 and mu1
    data_KbyS1 = bsxfun(@minus, data_KbyS1, mu1);
    data_KbyS1 = bsxfun(@rdivide, data_KbyS1, SD1);
    if(~isempty(data_KbyS2))
        if(size(data_KbyS1, 1) ~= size(data_KbyS2, 1))
            error('feature number for data_KbyS1 and data_KbyS2 must be the same!')
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
    index = logical(index_S1 + index_S2);
    data_KbyS = data_KbyS(:, index);
end

%data_KbyS = zscore(data_KbyS,0,2);
if(strcmp(type, 'corr'))
    FSM = CBIG_self_corr(data_KbyS);
elseif(strcmp(type,'Exponential'))
    FSM = exp(-1*scale*squareform(pdist(data_KbyS'))/K);
elseif(strcmp(type,'Gaussian'))
    FSM = exp(-1*scale*squareform(pdist(data_KbyS').^2)/K);
else
    FSM = CBIG_self_corr(data_KbyS);
end