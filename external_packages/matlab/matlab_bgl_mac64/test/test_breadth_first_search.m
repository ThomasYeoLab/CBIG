function distmap = test_breadth_first_search(A,u)

distmap = ipdouble(zeros(size(A,1),1));
distmap(u) = 0;
    
    function on_tree_edge(ei,u,v)
        distmap(v) = distmap(u)+1;
    end

breadth_first_search(A,u,struct('tree_edge',@on_tree_edge));

distmap = double(distmap);
end