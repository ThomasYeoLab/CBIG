function rval=rtest_1()

n = 49;
[A,b] = testmat(n,2);
x0 = [1:n]'/(n+1);
y0 = [1:n]'/(n+1);
x = repmat(x0,1,n);
y = repmat(y0',n,1);
xy = [x(:),y(:)];

rval = 0;
try    
    T = mst(A);
    T = T + diag(diag(A));
    rval = 1;
catch
    lasterr
end;

try
    A(1,2)= -1; 
    A(2,1)= -1; 
    T = prim_mst(A);
    rval = 0;
catch
    rval = 1;
end;

function [A,b] = testmat( n,stencil )
    h = 1/(n+1);
    
    % initialization
    x = [1:n]'/(n+1);
    y = [1:n]'/(n+1);
    u = zeros(n,n);
    
    % exact solution
    u0 = repmat(x.^4,1,n) + repmat(12*y'.^2,n,1);
    
    % matices
    T0 = -ones(n);
    T0 = sparse( triu(tril(T0,-1),-1) + triu(tril(T0,1),1) );
    I  = speye(n);
    if stencil == 1
        T = (-8*I - T0) / 3;
        B = (I - T0) / 3;
    else
        T = (-20*I - 4*T0) / 6;
        B = (4*I - T0) / 6;
    end
    
    A = kron(I,T) - kron(T0,B);

    % boundary
    b = zeros(n);
    if stencil == 1
        b(1,:) = b(1,:) + 12*y'.^2 + 12*(y-h)'.^2 + 12*(y+h)'.^2;
        b(n,:) = b(n,:) + 3 + 12*y'.^2 + 12*(y-h)'.^2 + 12*(y+h)'.^2;
        b(:,1) = b(:,1) + x.^2 + (x-h).^2 + (x+h).^2;
        b(:,n) = b(:,n) + x.^2 + (x-h).^2 + (x+h).^2 + 12*3;
    else
        b(1,:) = b(1,:) + 4* 12*y'.^2 + 12*(y-h)'.^2 + 12*(y+h)'.^2;
        b(n,:) = b(n,:) + 6 + 4* 12*y'.^2 + 12*(y-h)'.^2 + 12*(y+h)'.^2;
        b(:,1) = b(:,1) + 4* x.^2 + (x-h).^2 + (x+h).^2;
        b(:,n) = b(:,n) + 4* x.^2 + (x-h).^2 + (x+h).^2 + 12*6;
    end
    % fix corners
    b(1,n) = b(1,n) - 12;
    b(n,1) = b(n,1) - 1;
    b(n,n) = b(n,n) - 13;
    % normalize
    if stencil == 1
        b = - b / 3;
    else
        b = - b / 6;
    end
    % add source
    f = 12*x.^2 + 24;
    f = repmat(f,1,n);
    b = b + f * h^2;
    
    b = b(:);
return    