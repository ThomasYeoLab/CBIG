function M=caret_load(filename)
% M=caret_load(filename)
% loads a caret data file into a specific strcuture
% Types of possible files so far implemented: 
% .metric
% .surface_shape
% .paint
% .coord
% .topo
% .borderproj
% .latlon
%---------------------------------------------------
% v.1.0 Joern Diedrichsen 12/04 
% jdiedric@bme.jhu.edu
if(nargin<1)
    [file,path]=uigetfile('*.*','open file');
    filename=[path file];
end;
fid=fopen(filename,'r','ieee-be');
if (fid==-1)
    error(['Error: Could not find ' filename]);
    M=[];
    return;    
end;

% Figure out the type of file we are loading: 
s=strfind(filename,'.');
type=filename(s(end)+1:end);

% --------------------------------------------------------------------
% Load header of file 
% and analyze the header 
done_header=0;
i=1;
while (done_header==0)
    M.header{i}=fgetl(fid);
    if strcmp(M.header{i},'EndHeader')
        done_header=1;
    end;
    i=i+1;
end;
s=strmatch('encoding',char(M.header));
[dummy,M.encoding]=strread(M.header{s},'%s%s');

% --------------------------------------------------------------------
% Now read the data and format it, depending on the type of file 
switch (type)
    
    % --------------------------------------------------------------------
    % Coordinate file
    case {'coord'}
        if (strcmp(M.encoding,'ASCII'))
            M.data=textread(filename,'%f','headerlines',i-1,'delimiter',' ');    
            M.num_nodes=M.data(1);
            M.data=reshape(M.data(2:end),4,M.num_nodes)';
            M.index=M.data(:,1);
            M.data=M.data(:,2:end);
        else 
            M.num_nodes=fread(fid,1,'int32');
            M.data=fread(fid,inf,'float32');
            M.data=reshape(M.data,3,M.num_nodes)';
            M.index=[0:M.num_nodes-1];
        end;
    % --------------------------------------------------------------------
    % Coordinate file
    case {'latlon'}
        [M,t]=load_tags(fid,M);
        s=strmatch('tag-number-of-nodes',char(M.tags));
        [dummy,M.num_nodes]=strread(M.tags{s},'%s%d');
        if (strcmp(M.encoding,'ASCII'))
            M.data=textread(filename,'%f','headerlines',i-1,'delimiter',' ');    
            M.num_nodes=M.data(1);
            M.data=reshape(M.data(2:end),4,M.num_nodes)';
            M.index=M.data(:,1);
            M.data=M.data(:,2:end);
        else 
          %  M.num_nodes=fread(fid,1,'int64');
            M.data=fread(fid,inf,'float');
            M.data=reshape(M.data,4,M.num_nodes)';
            M.index=[0:M.num_nodes-1];
        end;
    
     % --------------------------------------------------------------------
    % topology files
    case 'topo'
        string=fgetl(fid);
        if (strcmp(string,'tag-version 1'))
            M.tag_version=1;
        else 
            M.tag_version=0;i=i-1;
        end;            
        if (strcmp(M.encoding,'ASCII'))
            buffer=textread(filename,'%d','headerlines',i,'delimiter',' '); 
            where=1;
            if (M.tag_version==0)
                M.num_nodes=buffer(where);where=where+1;
                for i=1:M.num_nodes
                    info=buffer(where:where+5);
                    M.Neighbor{info(1)+1}=buffer(where+7:2:where+info(2)*2+5)+1;
                    where=where+info(2)*2+6;
                    if(mod(i,1000)==0)
                        fprintf('.',i);
                    end;
                end;
            end;
            M.num_tiles=buffer(where);
            M.data=reshape(buffer(where+1:end),3,M.num_tiles)';
        else 
            M.data=fread(fid,inf,'int32');
            M.num_tiles=M.data(1);
            M.data=reshape(M.data(2:end),3,M.num_tiles)';
        end;
        M.data=M.data+1;  % Make Indices start at 1, not 0 
        
    % --------------------------------------------------------------------
    % metric files 
    case {'surface_shape','metric'}
        [M,t]=load_tags(fid,M);
        s=strmatch('tag-number-of-columns',char(M.tags));
        [dummy,M.num_cols]=strread(M.tags{s},'%s%d');
        s=strmatch('tag-number-of-nodes',char(M.tags));
        [dummy,M.num_rows]=strread(M.tags{s},'%s%d');
        s=strmatch('tag-column-name',char(M.tags));
        for k=1:length(s)
            S=strread(M.tags{s(k)},'%s',-1);
            name=horzcat(S{3:end});
            j=str2num(S{2});
            M.column_name{j+1}=name;
        end;
        s=strmatch('tag-column-color-mapping',char(M.tags));
        for k=1:length(s);
            [dummy,num,low,high]=strread(M.tags{s(k)},'%s %d %f %f',1);
            M.column_color_mapping(num+1,:)=[low high];
        end;
        if (strcmp(M.encoding,'ASCII'))
            M.data=textread(filename,'%f','headerlines',i+t-2,'delimiter',' ');    
            M.data=reshape(M.data,M.num_cols+1,M.num_rows)';
            M.index=M.data(:,1);
            M.data=M.data(:,2:end);
        else 
            M.data=fread(fid,inf,'float32');
            M.index=[1:M.num_rows]';
            M.data=reshape(M.data,M.num_cols,M.num_rows)';
        end;
    
    % --------------------------------------------------------------------
    % paint files     
    case 'paint'
        [M,t]=load_tags(fid,M);
        s=strmatch('tag-number-of-columns',char(M.tags));
        [dummy,M.num_cols]=strread(M.tags{s},'%s%d');
        s=strmatch('tag-number-of-nodes',char(M.tags));
        [dummy,M.num_rows]=strread(M.tags{s},'%s%d');
        s=strmatch('tag-number-of-paint-names',char(M.tags));
        [dummy,M.num_paintnames]=strread(M.tags{s},'%s%d');
        s=strmatch('tag-column-name',char(M.tags));
        for k=1:length(s)
            S=strread(M.tags{s(k)},'%s',-1);
            name=horzcat(S{3:end});
            j=str2num(S{2});
            M.column_name{j+1}=name;
        end;
        s=strmatch('tag-column-color-mapping',char(M.tags));
        for k=1:length(s);
            [dummy,num,low,high]=strread(M.tags{s(k)},'%s %d %f %f',1);
            M.column_color_mapping(num+1,:)=[low high];
        end;
        
        % read number of paintnames
        for pn=1:M.num_paintnames
            s=fgetl(fid);
            [num,M.paintnames(pn)]=strread(s,'%d %s');
        end;
        
        if (strcmp(M.encoding,'ASCII'))
            M.data=textread(filename,'%f','headerlines',i+t+pn-2,'delimiter',' ');    
            M.data=reshape(M.data,M.num_cols+1,M.num_rows)';
            M.index=M.data(:,1);
            M.data=M.data(:,2:end);
        else 
            M.data=fread(fid,inf,'int32');
            M.index=[1:M.num_rows]';
            M.data=reshape(M.data,M.num_cols,M.num_rows)';
        end;
        
    % --------------------------------------------------------------------
    % Border-projection files
    case 'borderproj'
        [M.numborders]=fscanf(fid,'%d',1);
        for b=1:M.numborders
            B.numborder=fscanf(fid,'%d',1);
            B.numpoints=fscanf(fid,'%d',1);
            B.name=fscanf(fid,'%s',1);
            B.p=fscanf(fid,'%f',4);
            B.c=fscanf(fid,'%f',3);
            data=fscanf(fid,'%f',[7 B.numpoints]);
            data=data';
            B.vertex=data(:,1:3);
            B.weight=data(:,5:7);
            M.Border(b)=B;
        end;
    otherwise 
        error(['caret_load: Unknown ot not implemented file type: ' type]);
end;
fclose(fid);


function [M,i]=load_tags(fid,M)
done_tags=0;
i=1;
while (done_tags==0)
    M.tags{i}=fgetl(fid);
    if strcmp(M.tags{i},'tag-BEGIN-DATA')
        done_tags=1;
    end;
    i=i+1;
end;
