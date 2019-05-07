function output_mtx = CBIG_MMLDA_cellstr2matrix(input_cell)

% output_mtx = CBIG_MMLDA_cellstr2matrix(input_cell)
% 
% Tranform a cell array to matrix.
%
% Input:
%   - input_cell    : input cell string
% 
% Output:
%   - output_mtx    : ouput matrix
%
% Example:
% input_cell = {'1', '2'};
% output_mtx = CBIG_MMLDA_cellstr2matrix(input_cell)
% output_mtx = [1, 2];
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

output_mtx = cell2mat(cellfun(@str2num, input_cell, 'uniformOutput', false));