function ind = CBIG_MMLDA_find_cell_in_cell(source_cell, target_cell)
% ind = CBIG_MMLDA_find_cell_in_cell(source_cell, target_cell)
%
% Find source cell array in target cell array and return the index in order.
% 
% Input:
%   - source_cell   : source cell array
%   - target_cell   : target cell array
%
% Output:
%   - ind           : index array 
%
% Example:
%   ind = CBIG_MMLDA_find_cell_in_cell({'abc'}, {'aaa', 'abc'});
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ind = zeros(1, length(source_cell));
for i = 1:length(source_cell)
    if sum(strcmp(target_cell, source_cell{i})) == 0
        ind(i) = 0;
    else
        ind(i) = find(strcmp(target_cell, source_cell{i}) == 1);
    end
end