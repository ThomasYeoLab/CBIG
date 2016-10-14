function caret_savemetric(filename,M)
% function caret_savemetric(filename,M)
% Save a M-file structure as a binary or ascii metric file 
% M.encoding should be {'ASCII','BINARY'}
% --------------------------------------------------------------
% v.1.0 Joern Diedrichsen 04/11/28
if nargin==1
    M=filename;
    filename='';
end;
if isempty(filename)
   [F,P]=uiputfile('*.metric','Save Metric file as');
   filename = [P,F];
end
fid=fopen(filename,'w','ieee-be');
if (fid==-1)
    fprintf('Error opening file %s\n',filename);
    return;
end;
linefeed=sprintf('\n');

% --------------------------------------------------------------
% Make Header: 
M.header={};
M.header{1}='BeginHeader';
M.header{2}=['encoding ' M.encoding{1}];
M.header{3}='EndHeader';
M.header{4}='tag-version 2';
M.header{5}=sprintf('tag-number-of-nodes %d',size(M.data,1));
M.header{6}=sprintf('tag-number-of-columns %d',size(M.data,2));
M.header{7}=sprintf('tag-title');
where=7;
for i=1:length(M.column_name)
    M.header{where+i}=sprintf('tag-column-name %d %s',i-1,M.column_name{i});
end;
where=where+length(M.column_name);
for i=1:length(M.column_name)
    M.header{where+i}=sprintf('tag-column-color-mapping %d %8.8f %8.8f',i-1,M.column_color_mapping(i,1),M.column_color_mapping(i,2));
end;
where=where+length(M.column_name);
M.header{where+1}='tag-BEGIN-DATA';

for line=1:length(M.header) 
    fprintf(fid,'%s\n',M.header{line});
end;
% Format string: 
if (strcmp(M.encoding,'ASCII'))
    M.data=[M.index M.data];
    format_str=['%i' repmat([' %6.6f'],1,size(M.data,2)-1) '\n'];
    fprintf(fid,format_str,M.data');
else
    fwrite(fid,M.data','float32');
end;
fclose(fid);
