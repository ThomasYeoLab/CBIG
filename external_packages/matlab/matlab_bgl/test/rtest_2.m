function rval=rtest_2()

rval = 0;
try
    load ../graphs/cs-stanford.mat
    [ci s] = components(A);
    if (max(ci) ~= 4391)
        error('Incorrect number of components');
    end;
    rval = 1;
catch
    lasterr
end