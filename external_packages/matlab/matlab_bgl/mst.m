function [out1 out2 out3] = mst(A,varargin)
% MST Compute a minimum spanning tree for an undirected graph A.
%
% There are two ways to call MST.
% T = mst(A)
% [i j v] = mst(A) 
% The first call returns the minimum spanning tree T of A.  
% The second call returns the set of edges in the minimum spanning tree.  
% The calls are related by 
%    T = sparse(i,j,v,size(A,1), size(A,1)); 
%    T = T + T';
% The optional algname parameter chooses which algorithm to use to compute
% the minimum spanning tree.  Note that the set of edges returned is not
% symmetric and the final graph must be explicitly symmetrized.
%
% This method works on undirected graphs.
%
% ... = mst(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.algname: the minimum spanning tree algorithm
%       ['prim' | {'kruskal'}]
%   options.edge_weight: a double array over the edges with an edge
%       weight for each node, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly, see Note 1.
%       [{'matrix'} | length(nnz(A)) double vector]
%   options.root: specify the root or starting vertex for the algorithm
%       This option only applies to prim's algorithm. 
%       [{'none'} | any vertex number]
%   options.fix_diag: remove any diagonal entries to get correct output
%       from Prim's algorithm [0 | {1}]; beware this option with the
%       edge_weight option too.
%
% Note: the input to this function must be symmetric, so this function
% ignores the 'notrans' default option and never transposes the input.
%
% Note 1: see EXAMPLES/REWEIGHTED_GRAPHS for how to reweight a symmetric
% graph correctly.  There are a few complicated details.
%
% Example:
%    load graphs/clr-24-1.mat
%    mst(A)
%
% See also PRIM_MST, KRUSKAL_MST

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-05-03: Changed to using kruskal as the default following problems
%    with prim due to negative edge weights.
%  2006-05-31: Added full2sparse option
%  2006-06-15: Fixed error with graph symmetric (T+T') instead of max(T,T')
%    found by Mark Cummins
%  2006-11-09: Temporary fix for disconnected graphs and the number of edges
%    in the mst is LESS than n-1.
%  2006-11-10: Added warning for prim with disconnected graphs.
%  2007-04-09: Fixed documentation typos.  (Thanks Chris Maes.)
%  2007-04-09: Fixed bug with 0 weighted graphs.  (Thanks Chris Maes.)
%  2007-04-20: Added edge weight option
%  2007-07-12: Fixed edge_weight documentation
%    Added note about symmetric edge weights
%  2007-12-14: Added rooted option for prim's algorithm
%  2008-10-07: Changed options parsing
%    Addressed issue with incorrect prim output and fixed matrix diagonal
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end
if trans, end % no trans check

options = struct('algname', 'kruskal', 'edge_weight', 'matrix', ...
    'root', 'none', 'fix_diag', 1);
options = merge_options(options,varargin{:});

fixed_diag= 0; 
if options.fix_diag && strcmp(options.algname,'prim') ... 
    && any(diag(A)), A = A - diag(diag(A)); fixed_diag= 1; end

% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix.
edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weights = 1;
    edge_weight_opt = options.edge_weight;
    if fixed_diag, warning('matlab_bgl:fix_diag',...
        'the diagonal was adjusted, the edge_weight option may be incorrect'); 
    end
end

if check
    % make sure the matrix is symmetric
    if ~edge_weights
        check_matlab_bgl(A,struct('sym',1,'values',1,...
            'noneg', strcmp(options.algname,'prim')));
    else
        check_matlab_bgl(A,struct());
        
        if strcmp(options.algname,'prim') && any(edge_weights < 0)
            error('matlab_bgl:invalidParameter', ...
                'the edge_weight array must be non-negative');
        end
        
        % check for symmetry
        [i j] = find(A);
        Av = sparse(i,j,edge_weight_opt,size(A,1), size(A,2));
        check_matlab_bgl(Av,struct('sym',1,'nodefault',1));

    end
    
    if strcmp(options.algname,'prim')
        if max(components(A)) > 1
            warning('mst:connected', ...
                ['The output from MST using Prim''s algorithm\n' ...
                 'on a disconnected graph is only a partial spanning tree.']);
        end
    end
    
end

% old temporary fix for incorrect number of edges
% num_components = max(components(A));

if strcmp(options.root,'none')
    root = 0; % a flag used to denote "no root" to the mex
elseif isa(options.root, 'double')
    root = options.root;
else
    error('matlab_bgl:invalidParameter', ...
        'options.root is not ''none'' or a vertex number.');
end

[i j v] = mst_mex(A,lower(options.algname),edge_weight_opt,root);

% old temporary fix for disconnected graphs
% if (num_components > 1)
%     i = i(1:end-(num_components-1));
%     j = j(1:end-(num_components-1));
%     v = v(1:end-(num_components-1));
% end

if (nargout == 1 || nargout == 0)
    T = sparse(i,j,v,size(A,1),size(A,1));
    T = T + T';
    out1 = T;
    
    if nnz(T) == 0 && nnz(A) > 0
        warning('mst:empty', ...
            ['MST is empty.  This can occur if you have reweighted\n' ...
             'the matrix with 0 edge weights.  Try the [i j v] = mst(...)\n' ...
             'call instead for that case.']);
    end
else
    out1 = i; 
    out2 = j;
    out3 = v;
end;





