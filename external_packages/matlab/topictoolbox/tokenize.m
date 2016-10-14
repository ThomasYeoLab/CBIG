function [ Query , ni ] = tokenize( inputquery )

s = inputquery;
ni = 0;
while length( s ) > 0
    [token, rem] = strtok(s);
    ni = ni + 1;
    Query{ ni } = token;
    
    s = rem;
end