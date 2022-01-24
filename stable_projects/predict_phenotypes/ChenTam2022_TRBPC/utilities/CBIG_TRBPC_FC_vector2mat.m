function FC_mat = CBIG_TRBPC_FC_vector2mat(FC_vector)

% FC_mat = CBIG_TRBPC_FC_vector2mat(FC_vector)
%
% This function covert vectorized lower-triangular FC back to a 419*419
% matrix
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

num_ROIs = 419;
if ~isequal(size(FC_vector),[num_ROIs*(num_ROIs-1)/2,1])
    error('size of FC vector must be 87571*1');
end

FC_mat = zeros(num_ROIs,num_ROIs);
index = 0;
for j = 1:(num_ROIs-1)
    for i = (j+1):num_ROIs
        index = index + 1;
        FC_mat(i,j) = FC_vector(index);
    end
    FC_mat(i,i) = 0; % Diagonal set to 0
end
for i = 1:num_ROIs
    for j = (i+1):num_ROIs
        FC_mat(i,j) = FC_mat(j,i);
    end
end