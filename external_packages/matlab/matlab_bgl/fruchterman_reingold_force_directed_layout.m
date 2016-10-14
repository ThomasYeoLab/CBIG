function X = fruchterman_reingold_force_directed_layout(A,varargin)
% FRUCHTERMAN_REINGOLD_FORCE_DIRECTED_LAYOUT A force directed graph layout
% 
% Compute the layout for an unweighted, undirected graph.  
% See
% http://www.boost.org/doc/libs/1_36_0/libs/graph/doc/fruchterman_reingold.html
% for information about the layout function and the parameters
%
% The temperature of the layout begins at the value initial_temp and 
% decreases to 0 over the number of iterations given. 
%
% ... = fruchterman_reingold_force_directed_layout(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.iterations: number of layout iterations [int > 0 | {100}]
%   options.initial_temp: starting temperature [double > 0 | {10}]
%   options.force_pairs: computation of forces between pairs of vertices can
%     either be approximated on a grid or computed between all pairs 
%     [{'grid'} | 'all']
%   options.width: width of layout area [double | {num_vertices(G)}]
%   options.height: height of layout area [double | {num_vertices(G)}]
%   options.progressive: whether to start from an old layout
%     [{0} | position matrix X]
%
% Note: this function does not depend on the non-zero values of A, 
% but only uses the non-zero structure of A
%
% Example:
%   G = grid_graph(6,5);
%   X = fruchterman_reingold_force_directed_layout(G);
%   gplot(G,X);
%   X = fruchterman_reingold_force_directed_layout(G,'initial_temp',100);
%
% See also KAMADA_KAWAI_SPRING_LAYOUT, GURSOY_ATUN_LAYOUT, LAYOUT

% David F. Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-09-25: Initial coding
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

n = num_vertices(A);
options = struct('iterations',100,'initial_temp',10,'force_pairs','grid',...
    'width',n,'height',n,'progressive',0);
options = merge_options(options,varargin{:});

if check, check_matlab_bgl(A,struct('sym',1)); end

progressive_opt = [];
if ~isscalar(options.progressive), progressive_opt = options.progressive; end

force_pair_type = 0;
switch options.force_pairs
    case 'grid', force_pair_type = 1;
    case 'all', force_pair_type = 0;
    otherwise, error('matlab_bgl:invalidParameter',...
            'force_pair = %s isn''t a recognized option',options.force_pair);
end

X = fruchterman_reingold_mex(A,options.iterations,options.initial_temp,...
    force_pair_type, options.width, options.height, progressive_opt);
