function avg = CBIG_AverageStructures(output_file, var_name, varargin)

% Average each variable of a structure across all mat files.
% 
%   avg = CBIG_AverageStructures(output_file, var_name, varargin)
%   Input:
%       output_file : output file name
%       var_name    : structure name
%       varargin    : all the mat files
%   Output:
%       avg         : average structure
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


for i = 1:length(varargin)
  load(varargin{i}, var_name);
  eval(['tmp = ' var_name ';']);
  if(i == 1)
    avg = tmp;
  else
    field_names = fieldnames(tmp);
    for k = 1:length(field_names)
        eval(['avg.' field_names{k} ' = avg.' field_names{k} ' + tmp.' field_names{k} ';']);
    end
  end
end

field_names = fieldnames(avg);
for k = 1:length(field_names)
   eval(['avg.' field_names{k} ' = avg.' field_names{k} '/' num2str(length(varargin)) ';']);
end

eval([var_name ' = avg;']);
save(output_file, var_name);
exit
