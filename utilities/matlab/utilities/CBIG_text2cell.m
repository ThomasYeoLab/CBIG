function [cell_array, num_lines] = CBIG_text2cell(text_file)

% [cell_array, num_lines] = CBIG_text2cell(text_file)
%
% This function reads from "text_file" and return the content of each line
% to be an entry of "cell_array".
% 
% Inputs:
%   - text_file
%     The input text file.
% 
% Outputs:
%   - cell_array
%     The output cell structure, where each entry corresponds to a line in
%     "text_file".
% 
%   - num_lines
%     The total number of lines in "text_file". It should equal to the
%     length of "cell_array".
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


num_lines = 0;
fid = fopen(text_file);
while (~feof(fid))
    num_lines = num_lines + 1;
    cell_array{num_lines} = fgetl(fid);
end
fclose(fid);


end

