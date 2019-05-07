function output_cell = CBIG_MMLDA_matrix2cellstr(input_mtx)

% output_cell = CBIG_MMLDA_matrix2cellstr(input_mtx)
% 
% Tranform an array to a cell array.
%
% Input:
%   - input_mtx     : input matrix
% 
% Output:
%   - output_cell   : ouput cell array
%
% Example:
% input_mtx = [1, 2];
% output_cell = CBIG_MMLDA_matrix2cellstr(input_mtx);
% output_cell = {'1', '2'};
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

input_mtx_cell = num2cell(input_mtx);
output_cell = cellfun(@num2str, input_mtx_cell, 'uniformOutput', false);