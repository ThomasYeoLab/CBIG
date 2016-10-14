function [d dt pred] = bfs_in_mbgl(A,u)
% BFS_IN_MBGL Reimplement the BFS function with MatlabBGL visitors.
%
% [d dt pred] = bfs_in_mbgl(A,u) 
%
% See BFS
%
% Example:
%    load ../graphs/bfs_example.mat
%    d = bfs_in_mbgl(A,1)


ip_d = ipdouble(-ones(num_vertices(A),1));
ip_dt = ipdouble(-ones(num_vertices(A),1));
ip_pred = ipdouble(zeros(1,num_vertices(A)));

ip_time = ipdouble(1);

    function time_discover_vertex(u)
        ip_dt(u) = ip_time(1);
        ip_time(1) = ip_time(1) + 1;
    end

    function distance_tree_edge(ei,u,v)
        ip_d(v) = ip_d(u)+1;
    end

    function pred_tree_edge(ei,u,v)
        ip_pred(v) = u;
    end

vis_distance = struct('tree_edge', @distance_tree_edge);
vis_time = struct('discover_vertex', @time_discover_vertex);
vis_pred = struct('tree_edge', @pred_tree_edge);

vis = combine_visitors(vis_distance, vis_time, vis_pred);

ip_d(u) = 0;
ip_dt(u) = ip_time(1);
breadth_first_search(A,u,vis);

d = double(ip_d);
dt = double(ip_dt);
pred = double(ip_pred);

% end the function
end