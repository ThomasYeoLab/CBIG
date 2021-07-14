function corr_mat = CBIG_TRBPC_FC_vector_2_mat(FC_vector)
% corr_mat = CBIG_TRBPC_FC_vector_2_mat(FC_vector)
%
% This function takes a vector of length 87571*1 and transforms it into a
% a 419 x 419 matrix
%
% required input:
% - FC_vector: a vector of length 87571*1
%
% Written by Jianzhong Chen & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_ROIs = 419;
if ~isequal(size(FC_vector),[num_ROIs*(num_ROIs-1)/2,1])
    error('size of FC vector must be 87571*1');
end

corr_mat = zeros(num_ROIs,num_ROIs);
index = 0;
for j = 1:(num_ROIs-1)
    for i = (j+1):num_ROIs
        index = index + 1;
        corr_mat(i,j) = FC_vector(index);
    end
    corr_mat(i,i) = 0; % Diagonal set to 0
end
for i = 1:num_ROIs
    for j = (i+1):num_ROIs
        corr_mat(i,j) = corr_mat(j,i);
    end
end