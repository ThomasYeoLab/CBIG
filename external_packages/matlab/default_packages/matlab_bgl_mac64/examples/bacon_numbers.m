function bn = bacon_numbers(A,u)
% BACON_NUMBERS Compute the Bacon numbers for a graph.
%
% bn = bacon_numbers(A,u) computes the Bacon numbers for all nodes in the 
% graph assuming that Kevin Bacon is node u.

% allocate storage for the bacon numbers
% the ipdouble call allocates storage that can be modified in place.
bn_inplace = ipdouble(zeros(num_vertices(A),1));

% implement a nested function that can refer to variables we declare.  In
% this case, we refer to the bn_inplace variable.  
function tree_edge(ei,u,v)
    bn_inplace(v) = bn_inplace(u)+1;
end

% setup the bacon_recorder visitor
bacon_recorder = struct();
bacon_recorder.tree_edge = @tree_edge;

% call breadth_first_search
breadth_first_search(A,u,bacon_recorder);

% convert the inplace storage back to standard Matlab storage to return.
bn = double(bn_inplace);

% the end line is required with nested functions to terminate the file
end