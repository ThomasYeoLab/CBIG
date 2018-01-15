function distmap = test_depth_first_search(A,u)

distmap = -1*ones(size(A,1),1);

distmap = ipdouble(distmap);
distmap(u) = 0;
    
    function on_tree_edge(ei,u,v)
        distmap(v) = distmap(u)+1;
    end

depth_first_search(A,u,struct('tree_edge',@on_tree_edge));

distmap = double(distmap);
end