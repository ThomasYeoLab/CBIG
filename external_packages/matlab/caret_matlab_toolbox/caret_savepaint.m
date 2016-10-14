function caret_savemetric(filename,M)
% function caret_savemetric(filename,M)
% save a paint file structure (as gotten from caret_load) 
% As a paint file readable by Caret (v 5.12)
if nargin==1
    M=filename;
    filename='';
end;
if isempty(filename)
   [F,P]=uiputfile('*.*','Save Cell Array as');
   filename = [P,F];
end
fid=fopen(filename,'w','ieee-be');
if (fid==-1)
    fprintf('Error opening file %s\n',filename);
    return;
end;
linefeed=sprintf('\n');
% Make Header: 
s=strmatch('encoding',char(M.header));
M.header{s}=sprintf('encoding %s',char(M.encoding));
M.tags={};
M.tags{1}='tag-version 1';
M.tags{2}=sprintf('tag-number-of-nodes %d',size(M.data,1));
M.tags{3}=sprintf('tag-number-of-columns %d',size(M.data,2));
M.tags{4}=sprintf('tag-number-of-paint-names %d',M.num_paintnames);
M.tags{5}=sprintf('tag-title');
where=5;
for i=1:length(M.column_name)
    M.tags{where+i}=sprintf('tag-column-name %d %s',i-1,M.column_name{i});
end;
where=where+length(M.column_name);
M.tags{where+1}='tag-BEGIN-DATA';

for line=1:length(M.header) 
    fprintf(fid,'%s\n',M.header{line});
end;
for line=1:length(M.tags) 
    fprintf(fid,'%s\n',M.tags{line});
end;

% Format string: 
for i=1:M.num_paintnames
    fprintf(fid,'%d %s\n',i-1,M.paintnames{i});
end;
if (strcmp(M.encoding,'ASCII'))
    M.data=[M.index M.data];
    format_str=['%i' repmat([' %6.6f'],1,size(M.data,2)-1) '\n'];
    fprintf(fid,format_str,M.data');
else
    fwrite(fid,M.data','int32');
end;
fclose(fid);
