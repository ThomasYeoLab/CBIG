function rval = rtest_6()

rval = 0;

try

    % create a line graph
    n = 10;
    A = sparse(1:n-1,2:n,1,n,n);
    A = A+A';

    u = 1; 
    v = 5;

    d = dist_uv(A,u,v);

    if any(d(v+1:end) > 0)
        error('breadth_first_search did not stop correctly');
    end
    
    rval = 1;
catch
    lasterr
end

end

function dmap = dist_uv(A,u,v)
    vstar = v;
    dmap = ipdouble(zeros(size(A,1),1));
    function stop=on_tree_edge(ei,u,v)
        dmap(v) = dmap(u)+1;
        stop = (v ~= vstar);
    end
    breadth_first_search(A,u,struct('tree_edge',@on_tree_edge));
    dmap = double(dmap);
end

