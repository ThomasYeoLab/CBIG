function CBIG_AverageCorrelationMatrixGeneral(output_file, varargin_text)

% Read a file list, average all mat files in the list and save output.
% 
%   CBIG_AverageCorrelationMatrixGeneral(output_file, varargin_text)
%   Input: 
%       output_file  : output file name
%       varargin_text: input list name
%   
%   Example:
%   clear
%   corr_mat = rand(100,100);
%   save ~/storage/test/1.mat corr_mat
%   corr_mat = rand(100,100);
%   save ~/storage/test/2.mat corr_mat
%   fid = fopen('~/storage/test/test.txt','w');
%   fprintf(fid, '~/storage/test/1.mat \n');
%   fprintf(fid, '~/storage/test/2.mat \n');
%   fclose(fid);
%   CBIG_AverageCorrelationMatrixGeneral('average', '~/storage/test/test.txt')
% 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% read in file lists
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
varargin

for i = 1:length(varargin)
    disp(varargin{i});
    x = load(varargin{i});
    
   if(i == 1)
       corr_mat = x.corr_mat;
   else
       corr_mat = corr_mat + x.corr_mat;
   end
end
corr_mat = corr_mat/length(varargin);
save(output_file, 'corr_mat');
