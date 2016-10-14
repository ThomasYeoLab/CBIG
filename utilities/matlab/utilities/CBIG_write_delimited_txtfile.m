function CBIG_write_delimited_txtfile(filename, headers_cell, alldata_str, out_info)

% Output specific field of headers cell array, all data string to textfile.
% 
%   CBIG_write_delimited_txtfile(filename, headers_cell, alldata_str, out_info)
%   Input:
%       filename    : filename
%       headers_cell: header cell array
%       alldata_str : all data string
%       out_info    : output field
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


fid = fopen(filename, 'w');

if(nargin > 3)
    index = [];
    for i = 1:length(headers_cell)
        if(sum(ismember(out_info, headers_cell{i})) > 0)
            index = [index; i];
        end
    end
    headers_cell = headers_cell(index);
    alldata_str = alldata_str(index);
end

% first write headers 
for i = 1:length(headers_cell)-1 
   fprintf(fid, '%s,', headers_cell{i}); 
end
fprintf(fid, '%s\n', headers_cell{end});

num_lines = length(alldata_str{1});
for i = 1:num_lines
   for j = 1:(length(headers_cell)-1)
       fprintf(fid, '%s,', alldata_str{j}{i});
   end
   fprintf(fid, '%s\n', alldata_str{end}{i});
end

fclose(fid);

