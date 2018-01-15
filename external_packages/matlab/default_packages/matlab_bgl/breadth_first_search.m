function breadth_first_search(A,u,bfs_visitor,varargin)
% BREADTH_FIRST_SEARCH Fully wrap the Boost breadth_first_search call
% including the bfs_visitor.
%
% breadth_first_search(A,u,vis) performs a breadth first traversal 
% of A starting from vertex u.  For each event defined by the bfs_visitor 
% structure below, the visitor is called with the either the name of the 
% vertex (u), or the edge index and it's source and target (ei,u,v).  
%
% See http://www.boost.org/libs/graph/doc/BFSVisitor.html for a description
% of the events.
% 
% bfs_visitor is a struct with the following optional fields
%    vis.initialize_vertex(u)
%    vis.discover_vertex(u)
%    vis.examine_vertex(u)
%    vis.examine_edge(ei,u,v)
%    vis.tree_edge(ei,u,v)
%    vis.non_tree_edge(ei,u,v)
%    vis.gray_target(ei,u,v)
%    vis.black_target(ei,u,v)
%    vis.finish_vertex(u)
% Each bfs_visitor parameter should be a function pointer, which returns 0
% if the bfs should stop.  (If the function does not return anything, the
% bfs continues.)
%
% This method works on directed graphs.
% The runtime is O(V+E), excluding the complexity of the visitor
% operations.
%
% Realistically, this function must be used with the
% pass-by-reference/in-place modification library.  
%
% ... = breadth_first_search(A,u,vis,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   There are no additional options for this function.
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
%     end;
%     breadth_first_search(A,u,struct('tree_edge',@on_tree_edge));
%   end;
%
% See also BFS

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-05-15: Initial version
%  2006-05-31: Added full2sparse check 
%  2007-04-17: Fixed documentation
%% 

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

if check
    % no additional input checks
    check_matlab_bgl(A,struct());
end

if trans
    A = A';
end

%bfs_visitor = merge_structs(bfs_visitor, empty_bfs_visitor);

% The 101 is the flag for calling BFS, not DFS
bfs_dfs_vis_mex(A,u,bfs_visitor,101);


