function out = CBIG_TRBPC_load_mat(mat_file)

% out = CBIG_TRBPC_load_mat(mat_file)
%
% This function load a variable from a .mat file without knowing the variable name
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

mat_struct = load(mat_file);
names = fieldnames(mat_struct);
if length(names) ~= 1
    erorr('one and only one variable expectexd in the .mat file')
end
varname = names{1};
out = mat_struct.(varname);

end