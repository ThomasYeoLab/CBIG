function rval=rtest_7_karsi
% Test the flow-graph that Karci identified as causing a problem with the
% min-cut computation.  The min-cut value should always be equal to the
% max-flow value.  
G =  [0    46     0    46     0     0     0     0     0     0   461;
    0     0    46     0    46     0     0     0     0     0   236;
    0     0     0     0     0    46     0     0     0     0    23;
    0     0     0     0    46     0    46     0     0     0   300;
    0     0     0     0     0    46     0    46     0     0    46;
    0     0     0     0     0     0     0     0    46     0    46;
    0     0     0     0     0     0     0    46     0     0    23;
    0     0     0     0     0     0     0     0    46     0    73;
    0     0     0     0     0     0     0     0     0     0   183;
    46    46    56    46    95    58    56    23     0     0     0;
    0     0     0     0     0     0     0     0     0     0     0];
G=sparse(G);
try
    [f c] = max_flow(G,10,11);
    [i j v] = find(G);
    I = c(i)>c(j); % form an indicator variable for edges crossing the cut
    cv = sum(I.*v);
    if cv~=f, error('Incorrect cut value'); end
    rval=1;
catch 
    lasterr
end
