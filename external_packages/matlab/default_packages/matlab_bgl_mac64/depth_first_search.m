function depth_first_search(A,u,dfs_visitor,varargin)
% DEPTH_FIRST_SEARCH Fully wrap the Boost depth_first_search call
% including the dfs_visitor.
%
% depth_first_search(A,u,dfs_visitor) performs a depth first traversal 
% of A starting from vertex u.  For each event defined by the dfs_visitor 
% structure below, the visitor is called with the either the name of the 
% vertex (u), or the edge index and it's source and target (ei,u,v).  
%
% See http://www.boost.org/libs/graph/doc/DFSVisitor.html for a description
% of the events.
% 
% dfs_visitor is a struct with the following optional fields
%    vis.initialize_vertex(u)
%    vis.start_vertex(u)
%    vis.discover_vertex(u)
%    vis.examine_edge(ei,u,v)
%    vis.tree_edge(ei,u,v)
%    vis.back_edge(ei,u,v)
%    vis.forward_or_cross_edge(ei,u,v)
%    vis.finish_vertex(u)
% Each dfs_visitor parameter should be a function pointer, which returns 0
% if the dfs should stop.  (If the function does not return anything, the
% dfs continues.)
%
% This method works on directed graphs.
% The runtime is O(V+E), excluding the complexity of the visitor
% operations.
%
% Realistically, this function must be used with the
% pass-by-reference/in-place modification library.  
%
% ... = depth_first_search(A,u,vis,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.full: compute the full dfs instead of the dfs of
%      the current component (see Note 1) [{0} | 1]
%
% Note 1: When computing the full dfs, the vertex u is ignored, vertex 1 is
% always used as the starting vertex.  
%
% Note: this function does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.
%
% Example:
%   This example finds the distance to a single point and stops the search.
%   function dist_uv(A,u,v)
%     vstar = v;
%     dmap = ipdouble(zeros(size(A,1),1));
%     function stop=on_tree_edge(ei,u,v)
%       dmap(v) = dmap(u)+1;
%       stop = (v ~= vstar);
%     end
%     depth_first_search(A,u,struct('tree_edge',@on_tree_edge));
%   end
%
% See also DFS

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2006-05-21: Initial version
%  2006-05-31: Added full2sparse check
%  2007-07-24: Fixed example
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

if check
    % no additional input checks
    check_matlab_bgl(A,struct());
end

if trans, A = A'; end

% parse the optional parameters
full = 0;
if ~isempty(varargin)
    optionsu = merge_options(struct(),varargin{:});
    if (isfield(optionsu,'full'))
        full = optionsu.full;
    end
end

if full
    % 202 is the call for dfs with full searches
    bfs_dfs_vis_mex(A,u,dfs_visitor,202);
else
    % 201 is the call for dfs with partial searches
    bfs_dfs_vis_mex(A,u,dfs_visitor,201);
end


