function CBIG_PartitionSubjectsVarargin(num_subjects, varargin_txt, out_path)

% Split subject's list into several lists and each list contains num_subjects.
% 
%   CBIG_PartitionSubjectsVarargin(num_subjects, varargin_txt, out_path)
%   Input:
%       num_subjects: number of subjects included in each text file
%       varargin_txt: input subject's list
%       out_path    : output path
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(num_subjects))
   num_subjects = str2num(num_subjects); 
end

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

start = 1;
stop = num_subjects;
count = 0;
while(stop < length(varargin))
    
    count = count + 1;
    fid = fopen([out_path '_' num2str(num_subjects) '_' num2str(count) '.txt'], 'w');
    for j = start:stop
        fprintf(fid, '%s\n', varargin{j});
    end
    fclose(fid);
    
    start = start + num_subjects;
    stop = stop + num_subjects;
end

exit