function caret_save(filename,M)
% function caret_save(filename,M)
% Save a M-file structure as a binary or ascii metric, coord, topo, pait file 
% M.encoding should be {'ASCII','BINARY'}
% --------------------------------------------------------------
% v.1.0 Joern Diedrichsen 04/11/28
if nargin==1
    M=filename;
    filename='';
end;
if isempty(filename)
    [F,P]=uiputfile('*.*','Save Structure as');
    filename = [P,F];
end
fid=fopen(filename,'w','ieee-be');
if (fid==-1)
    fprintf('Error opening file %s\n',filename);
    return;
end;
linefeed=sprintf('\n');

% Figure out the type of file we are loading: 
s=strfind(filename,'.');
type=filename(s(end)+1:end);

% --------------------------------------------------------------------
% Now switch the format and make the header depending on the type of file 
switch (type)
    % --------------------------------------------------------------------
    % Coordinate file
    case 'coord'
        if (~isfield(M,'header'))
            M.header{1}='BeginHeader';
            M.header{2}=['encoding ' M.encoding{1}];
            M.header{3}='EndHeader';
        end;
        for line=1:length(M.header) 
            fprintf(fid,'%s\n',M.header{line});
        end;
        % Format string: 
        if (strcmp(M.encoding,'ASCII'))
            fprintf(fid,'%d\n',M.num_nodes);
            M.data=[M.index M.data];
            format_str=['%i %3.3f %3.3f %3.3f\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.num_nodes,'int32');
            fwrite(fid,M.data','float32');
        end;
    % --------------------------------------------------------------------
    % Coordinate file
    case 'topo'
        if (~isfield(M,'header'))
        M.header{1}='BeginHeader';
        M.header{2}=['encoding ' M.encoding{1}];
        M.header{3}='filetype topo';
        M.header{4}='perimeter_id CUT';
        M.header{5}='EndHeader';
        end;
        for line=1:length(M.header) 
            fprintf(fid,'%s\n',M.header{line});
        end;
        M.data=M.data-1;
        if (strcmp(M.encoding,'ASCII'))
            format_str=['%d %d %d\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.num_tiles,'int32');
            fwrite(fid,M.data','int32');
        end;


    % --------------------------------------------------------------
    % Metric file 
    case 'metric'
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
end;
fclose(fid);
