ext_path = pwd;
addpath(ext_path);

path_to_add = fullfile(ext_path, 'minka');
if(exist(path_to_add, 'dir'))
   addpath(genpath(path_to_add)); 
end

path_to_add = fullfile(ext_path, 'topictoolbox');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_path, 'others');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

if(strcmp(mexext, 'mexmaci64'))
	path_to_add = fullfile(ext_path, 'matlab_bgl_mac64');
else
	path_to_add = fullfile(ext_path, 'matlab_bgl');
end
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_path, 'gifti_toolbox');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_path, 'edit_distances');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_path, 'caret_matlab_toolbox');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_path, 'DSP');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_path, 'cifti-matlab');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end
