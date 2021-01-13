function output = err(n,maxit,pmax,family)

if n==0
    output.n=0;
    output.fatal=false;
    output.msg='';
   
else
    switch family
        case 'gaussian'
            output = err_elnet(n,maxit,pmax);
        case 'binomial'
            output = err_lognet(n,maxit,pmax);
        case 'multinomial'
            output = err_lognet(n,maxit,pmax);
        case 'poisson'
            output = err_fishnet(n,maxit,pmax);
        case 'cox'
            output = err_coxnet(n,maxit,pmax);
        case 'mrelnet'
            output = err_mrelnet(n,maxit,pmax);
    end
    output.msg = sprintf('from glmnet Fortran code (error code %d); %s', n, output.msg);
end


%------------------------------------------------------------------
% End private function err
%------------------------------------------------------------------

function output = err_elnet(n,maxit,pmax)

if (n > 0)  %fatal error
    if (n < 7777)
        msg = 'Memory allocation error; contact package maintainer';
    elseif (n == 7777)
        msg = 'All used predictors have zero variance';
    elseif (n == 10000)
        msg = 'All penalty factors are <= 0';
    else
        msg = 'Unknown error';
    end
    output.n = n;
    output.fatal = true;
    output.msg = msg;
elseif (n < 0)  %non-fatal error
    if (n > -10000)
        msg = sprintf('Convergence for %dth lambda value not reached after maxit=%d iterations; solutions for larger lambdas returned',-n,maxit);
    end
    if (n < -10000)
        msg = sprintf('Number of nonzero coefficients along the path exceeds pmax=%d at %dth lambda value; solutions for larger lambdas returned',pmax,-n-10000);
    end
    output.n = n;
    output.fatal = false;
    output.msg = msg;
end


function output = err_lognet(n,maxit,pmax)

output = err_elnet(n,maxit,pmax);
if (n < -20000)
    output.msg = sprintf('Max(p(1-p),1.0e-6 at %dth value of lambda; solutions for larger values of lambda returned',-n-20000);
end
if ~strcmp(output.msg, 'Unknown error')
    return;
end
if (8000 < n) && (n < 9000)
    msg = sprintf('Null probability for class%d< 1.0e-5', n-8000);
elseif (9000 < n) && (n < 10000)
    msg = sprintf('Null probability for class%d> 1.0 - 1.0e-5',n-9000);
else
    msg = 'Unknown error';
end
output.n = n;
output.fatal = true;
output.msg = msg;


function output = err_fishnet(n,maxit,pmax)

output = err_elnet(n,maxit,pmax);
if ~strcmp(output.msg, 'Unknown error')
    return;
end
if (n == 8888)
    msg = 'Negative response values - should be counts';
elseif (n == 9999)
    msg = 'No positive observation weights';
else
    msg = 'Unknown error';
end
output.n = n;
output.fatal = true;
output.msg = msg;


function output = err_coxnet(n,maxit,pmax)

if (n > 0)  %fatal error
    output = err_elnet(n,maxit,pmax);
    if ~strcmp(msg, 'Unknown error')
        return;
    end
    if (n == 8888)
        msg = 'All observations censored - cannot proceed';
    elseif (n == 9999)
        msg = 'No positive observation weights';
    elseif (n == 20000) || (n == 30000)
        msg = 'Inititialization numerical error; probably too many censored observations';
    else
        msg = 'Unknown error';
    end
    output.n = n;
    output.fatal = true;
    output.msg = msg;
elseif (n < 0)
    if (n <= -30000)
        msg = sprintf('Numerical error at %dth lambda value; solutions for larger values of lambda returned',-n-30000);
        output.n = n;
        output.fatal = false;
        output.msg = msg;
    else
        output = err_elnet(n,maxit,pmax);
    end
end


function output = err_mrelnet(n,maxit,pmax)

if (n == 90000)
    msg = 'Newton stepping for bounded multivariate response did not converge';
    output.n = n;
    output.fatal = false;
    output.msg = msg;
else
    output = err_elnet(n,maxit,pmax);
end