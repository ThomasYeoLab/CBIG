function FC_vector = CBIG_FC_mat2vector(FC_mat)

% FC_vector = CBIG_FC_mat2vector(FC_mat)
%
% This function coverts a FC matrix into the vectorized lower triangle.
% Matlab extracts the lower triangle row by row. For example, for matrix
% [1, 0.5, 0.6, 0.7; 0.5, 1, -0.5, 0.2; 0.6, -0.5, 1 ,0.3; 0.7, 0.2, 0.3 ,1],
% passing in the lower triangle indices to convert it to vector form will
% extract [0.5, 0.6, 0.7], and then [-0.5, 0.2] and then [0.3] and concatenate
% them into a vector in that order.
%
% Inputs:
%
% - FC_mat
%   A #edges x #edges FC matrix
%
% Outputs:
%
% - FC_vector
%   A #edges x 1 vector of the edges from the FC matrix.
%
% Example (Converting a full 4x4 FC matrix to the lower triangle):
% FC_vector = CBIG_FC_mat2vector( ...
%                [1, 0.5, 0.6, 0.7; 0.5, 1, -0.5, 0.2; 0.6, -0.5, 1 ,0.3; 0.7, 0.2, 0.3 ,1])
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% check if dimensions are correct
FC_size = size(FC_mat);
if ~isequal(FC_size(1),FC_size(2))
    error('Must be a square matrix!');
end

% get indices to extract
FC_indices = ones(FC_size);
FC_indices = logical(tril(FC_indices, -1));
% extract lower triangle
FC_vector = FC_mat(FC_indices);

end