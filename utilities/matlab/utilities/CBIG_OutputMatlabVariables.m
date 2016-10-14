function CBIG_OutputMatlabVariables(matlab_data_file, eval_string, output_txt_file, prefix_string)

% Output a matlab variable to a text file.
% 
%   CBIG_OutputMatlabVariables(matlab_data_file, eval_string, output_txt_file, prefix_string)
%   Input: 
%       matlab_data_file    : a mat file
%       eval_string         : variable that you want to output
%       output_txt_file     : output text file
%       prefix_string       : prefix before the variable in text file
%   Example:
%   CBIG_OutputMatlabVariables('something.mat', 'metrics.nodal_pathlengths', 'tmp.txt', prefix_string) 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


load(matlab_data_file);

variable = eval(eval_string);
fid = fopen(output_txt_file, 'w');

if(nargin == 4)
    fprintf(fid, '%-40s = %-12.5f\n', prefix_string, variable);
else
    fprintf(fid, '%-12.5f\n', variable);
end
fclose(fid);

exit