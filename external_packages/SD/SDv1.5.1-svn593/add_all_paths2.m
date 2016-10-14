if(ispc)
    check_dir = 'dir';
else
    check_dir = 'ls';
end

curr_path = pwd;

addpath(pwd);

[status, y] = system([check_dir ' AnalysisTools']);
if(status == 0)
   addpath(fullfile(curr_path, 'AnalysisTools')); 
end

[status, y] = system([check_dir ' AtlasSharpnessAndLabelMatch']);
if(status == 0)
   addpath(fullfile(curr_path, 'AtlasSharpnessAndLabelMatch')); 
end

[status, y] = system([check_dir ' BasicTools']);
if(status == 0)
   addpath(fullfile(curr_path, 'BasicTools')); 
end

[status, y] = system([check_dir ' SphericalDemons']);
if(status == 0)
   addpath(fullfile(curr_path, 'SphericalDemons'));
   addpath(fullfile(curr_path, 'SphericalDemons', 'freesurfer'));
end

[status, y] = system([check_dir ' FeatureSelection']);
if(status == 0)
   addpath(fullfile(curr_path, 'FeatureSelection')); 
end

[status, y] = system([check_dir ' FuncCoord']);
if(status == 0)
   addpath(fullfile(curr_path, 'FuncCoord')); 
end

[status, y] = system([check_dir ' kd_tree']);
if(status == 0)
   addpath(fullfile(curr_path, 'kd_tree'));
end

[status, y] = system([check_dir ' MARS2']);
if(status == 0)
   addpath(fullfile(curr_path, 'MARS2'));
end

[status, y] = system([check_dir ' libsvm-mat-2.86-1']);
if(status == 0)
   addpath(fullfile(curr_path, 'libsvm-mat-2.86-1'));
end

[status, y] = system([check_dir ' overcomplete_wavelets']);
if(status == 0)
   addpath(fullfile(curr_path, 'overcomplete_wavelets'));
   addpath(fullfile(curr_path, 'overcomplete_wavelets', 'yawtb_interface'));
end
