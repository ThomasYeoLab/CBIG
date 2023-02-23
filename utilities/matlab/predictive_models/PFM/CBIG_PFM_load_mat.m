function out = CBIG_PFM_load_mat(mat_file)

% out = CBIG_PFM_load_mat(mat_file)
%
% This function load a variable from a .mat file without knowing the variable name
%
% Inputs:
%   - mat_file
%     Data path for mat file
% 
% Outputs:
%   - out
%     Content of mat file
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

mat_struct = load(mat_file);
names = fieldnames(mat_struct);
if length(names) ~= 1
    error('one and only one variable expected in the .mat file')
end
varname = names{1};
out = mat_struct.(varname);

end