
curr_path = pwd;
addpath(pwd);

if(exist('AnalysisTools', 'dir'))
	addpath(fullfile(curr_path, 'AnalysisTools'));
end

if(exist('AtlasSharpnessAndLabelMatch', 'dir'))
   addpath(fullfile(curr_path, 'AtlasSharpnessAndLabelMatch')); 
end

if(exist('BasicTools', 'dir'))
   addpath(fullfile(curr_path, 'BasicTools')); 
end

if(exist('SphericalDemons', 'dir'))
   addpath(fullfile(curr_path, 'SphericalDemons'));
   addpath(fullfile(curr_path, 'SphericalDemons', 'freesurfer'));
end

if(exist('FeatureSelection', 'dir'))
   addpath(fullfile(curr_path, 'FeatureSelection')); 
end

if(exist('FuncCoord', 'dir'))
   addpath(fullfile(curr_path, 'FuncCoord')); 
end

if(exist('kd_tree', 'dir'))
   addpath(fullfile(curr_path, 'kd_tree'));
end

if(exist('MARS2', 'dir'))
   addpath(fullfile(curr_path, 'MARS2'));
end

if(exist('libsvm-mat-2.86-1', 'dir'))
   addpath(fullfile(curr_path, 'libsvm-mat-2.86-1'));
end

if(exist('overcomplete_wavelets', 'dir'))
   addpath(fullfile(curr_path, 'overcomplete_wavelets'));
   addpath(fullfile(curr_path, 'overcomplete_wavelets', 'yawtb_interface'));
end
