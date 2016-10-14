function rval=rtest_5_henderson()

rval = 0;
try
    rand('state',0);
    A = erdos_reyni(100,.1);
    if full(max(max(A))) > 1
        error('erdos_reyni.m does not return an adjacency matrix.');
    end
    rval = 1;
catch
    lasterr
end