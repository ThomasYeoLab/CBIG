function CBIG_CreateRandomVarargin(num_subjects, num_randomness, varargin_txt, out_path)

% Create num_randomness list of files, each containing num_subjects names
% 
%   CBIG_CreateRandomVarargin(num_subjects, num_randomness, varargin_txt, out_path)
%   Input:
%       num_subjects    : number of subjects
%       num_randomness  : number of randomizations
%       varargin_txt    : input list
%       out_path        : output path
%   Output:
%       the output are some txt files: out_path '/rand_' num2str(num_subjects) '_' num2str(random) '.txt'
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(ischar(num_subjects))
   num_subjects = str2num(num_subjects); 
end

if(ischar(num_randomness))
   num_randomness = str2num(num_randomness); 
end

rand('twister',sum(100*clock))

% read in file lists
fid = fopen(varargin_txt, 'r');
i = 0;
while(1);
   tmp = fscanf(fid, '%s\n', 1);
   if(isempty(tmp))
       break
   else
       i = i + 1;
       varargin{i} = tmp;
   end
end
fclose(fid);

for i = 1:num_randomness
   
    index = randperm(length(varargin));
    fid = fopen([out_path '/rand_' num2str(num_subjects) '_' num2str(i) '.txt'], 'w');
    for j = 1:num_subjects
        fprintf(fid, '%s\n', varargin{index(j)});
    end
    fclose(fid);
end