function CBIG_AvgFreeSurferVolumes(varargin_text, output_file)

% Average all freesurfer volumes.
% 
%   CBIG_AvgFreeSurferVolumes(varargin_text, output_file)
%   Input:
%       varargin_text: text file including all mat files
%       output_file  : output file name
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% read in files
fid = fopen(varargin_text, 'r');
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

tic
for i = 1:length(varargin)
  disp([num2str(i) ': ' varargin{i}]);
    x = MRIread(varargin{i});
    if(sum(isnan(x.vol(:))) > 0)
        disp(['Warning: ' varargin{i} ' contains ' num2str(sum(isnan(x.vol(:)))) ' isnan .']);
    end
    
    if(i == 1)
        output = x;
    else
        output.vol = output.vol + x.vol;
    end
end
toc
clear x;
output.vol = output.vol/length(varargin);
MRIwrite(output, output_file);
%exit
