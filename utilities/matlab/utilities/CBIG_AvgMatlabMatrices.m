function CBIG_AvgMatlabMatrices(varargin_text, input_var_name, output_file, output_var_name, dimcount, exit_flag)

% Average one variable of all mat files in one text file.
%   
%   CBIG_AvgMatlabMatrices(varargin_text, input_var_name, output_file, output_var_name, dimcount)
%   Input:
%       varargin_text   : text file including all mat files
%       input_var_name  : variable name that you want to average
%       output_file     : output file name
%       output_var_name : output variable name
%       dimcount        : 0 = average matrix across all mat files; 1 = average each column of matrix across all mat files;
%                         2 = average each row of matrix across all mat files
%       exit_flag       : 0 = don't exit matlab when this function is done;
%                         1 = exit matlab when this function is done;
%       
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(nargin < 6)
   exit_flag = 1; 
else
   if(ischar(exit_flag))
      exit_flag = str2num(exit_flag); 
   end
end

if(nargin < 5)
   dimcount = 0; 
else
   if(ischar(dimcount))
      dimcount = str2num(dimcount); 
   end
end

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
    disp(num2str(i));
    x = load(varargin{i});
    eval(['x = x.' input_var_name ';']);
    if(sum(isnan(x(:))) > 0)
        disp(['Warning: ' varargin{i} ' contains ' num2str(sum(isnan(x(:)))) ' isnan .']);
    end

    if(i == 1)
        output = x;
        count = zeros(size(x));
    else
        output = output + x;
    end

    if(dimcount > 0)
        if(length(size(x)) > 2 || dimcount > 2)
            error('Does not handle matrices of dimensions higher than 2');
        end
        tmp_count = squeeze(sum(abs(x), dimcount));
        
        if(dimcount == 1)
            count = count + repmat(double(tmp_count~=0), size(x, 1), 1);
        else
            count = count + repmat(double(tmp_count~=0), 1, size(x, 2));
        end
    end
end
toc
clear x;
if(dimcount == 0)
    output = output/length(varargin);
else
    output = output./count;
end
eval([output_var_name ' = output;']);
save(output_file, output_var_name,'-v7.3');

%save('count.mat', 'count');
if (exit_flag == 1)
    exit
end