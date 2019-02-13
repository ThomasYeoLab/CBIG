function CBIG_cell2text(cell_var, filename)

% CBIG_cell2text(cell_var, filename)
% 
% This function writes a cell array (cell_var) into a text file "filename".
% 
% Inputs:
%   - cell_var
%     The cell array that the user wants to write out.
% 
%   - filename
%     Full path of the output text file.
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


fid = fopen(filename, 'w');
formatSpec = '%s\n';

for row = 1:length(cell_var);
    fprintf(fid, formatSpec, cell_var{row});
end

fclose(fid);