function FC_mat = CBIG_FC_vector2mat(FC_vector)

% FC_mat = CBIG_FC_vector2mat(FC_vector)
%
% This function coverts a vectorized lower-triangular FC back to a num_ROIs*num_ROIs
% matrix. num_ROIs is automatically calculated from the length of the FC_vector. The order 
% in which the matrix is filled corresponds to the way matlab extracts the lower 
% triangle (row by row) (See CBIG_FC_mat2vector for details).
%
% Inputs:
%
% - FC_vector
%   Vectorized lower-triangular FC (#edges x 1)
%
% Outputs:
%
% - FC_mat
%   A #edges x #edges FC matrix
%
% Example (Converting the lower triangle from a 4x4 FC matrix to the full matrix):
% FC_mat = CBIG_FC_mat2vector([0.5, 0.6, 0.7, -0.5, 0.2,0.3])
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% compute dimensions for final matrix
num_ROIs = (1 + sqrt(1 + (4*2*size(FC_vector,1))))/ 2;
% check whether num_ROIs is a whole number, print final size so user can
% verify if it is correct
if mod(num_ROIs,1) ~= 0
    error('Vector length incorrect: cannot be reshaped into square matrix');
end
disp(['Converting to '  num2str(num_ROIs) ' x ' num2str(num_ROIs) ' matrix'])

FC_mat = zeros(num_ROIs,num_ROIs);
index = 0;
% fill in lower triangle
for j = 1:(num_ROIs-1)
    for i = (j+1):num_ROIs
        index = index + 1;
        FC_mat(i,j) = FC_vector(index);
    end
    FC_mat(j,j) = 1; % Diagonal set to 1
end
FC_mat(num_ROIs,num_ROIs) = 1; % Diagonal set to 
% mirror for upper triangle
for i = 1:num_ROIs
    for j = (i+1):num_ROIs
        FC_mat(i,j) = FC_mat(j,i);
    end
end

end