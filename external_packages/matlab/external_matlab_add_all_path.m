ext_path = pwd;
ext_default_path = fullfile(ext_path, 'default_packages');

if(strcmp(mexext, 'mexmaci64'))
	path_to_add = fullfile(ext_default_path, 'matlab_bgl_mac64');
else
	path_to_add = fullfile(ext_default_path, 'matlab_bgl');
end
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'edit_distances');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'DSP');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'cifti-matlab');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'mtimesx_20110223');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'transforms');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'io');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'stats');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'others');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'figure_utilities');
if(exist(path_to_add, 'dir'))
   addpath(path_to_add); 
end

path_to_add = fullfile(ext_default_path, 'graph_cut');
if(exist(path_to_add, 'dir'))
   addpath(genpath(path_to_add)); 
end

path_to_add = fullfile(ext_default_path, 'WashU_gradients');
if(exist(path_to_add, 'dir'))
   addpath(genpath(path_to_add)); 
end

path_to_add = fullfile(ext_default_path, 'FDR');
if(exist(path_to_add, 'dir'))
   addpath(genpath(path_to_add)); 
end
