function CBIG_AverageRuns(output, varname, varargin)

% Average one variable (matrix) of all runs
% 
%   CBIG_AverageRuns(output, varname, varargin)
%   Input:
%       output  : name of output mat file
%       varname : variable name that you want to average 
%       varargin: all the mat files
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


for i = 1:length(varargin)
  load(varargin{i}, varname); % load varname variable in a mat file
  if(i == 1)
    eval(['tmp = ' varname '; ']);
  else
    eval(['tmp = tmp + ' varname '; ']);
  end
end
tmp = tmp/length(varargin);
data = tmp;
save(output, 'data');
exit