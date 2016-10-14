function M=caret_struct(kind,varargin);
% function M=caret_struct(kind,varargin);
%   Generates a new structure for caret-files 
%   Examples: 
%   M=caret_struct('coord','data',data); 
%   T=caret_struct('metric','data',data);
%   ...
c=1;
M.encoding={'BINARY'};
while c<=length(varargin)
    switch varargin{c}
        case {'data','num_nodes','num_cols','encoding','column_name'}
            M=setfield(M,varargin{c},varargin{c+1});
            c=c+2;
    end;
end;

switch (kind)
    case 'coord'
        if (isfield(M,'data'))
            M.num_nodes=size(M.data,1);
            if size(M.data,2)~=3
                error('Coord-struct needs three columns');
            end;
        else
            error('Coord-struct needs data');
        end;
    case 'metric'
        if (isfield(M,'data'))
            M.num_rows=size(M.data,1);
            M.num_cols=size(M.data,2);
            if (~isfield(M,'column_name'))
                for c=1:M.num_cols
                    M.column_name{c}=sprintf('Column_%2.2d',c);    
                end;
            end;
            if (~isfield(M,'column_color_mapping'))
                for c=1:M.num_cols
                    M.column_color_mapping(c,1:2)=[-5 5];
                end;
            end;
            
            
        else
            error('Metric-struct needs data');
        end;
end;