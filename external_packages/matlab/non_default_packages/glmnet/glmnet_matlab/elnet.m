function fit = elnet(x, is_sparse, irs, pcs, y, weights, offset, gtype, ...
    parm, lempty, nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, ...
    isd, intr, maxit, family)

ybar = y' * weights/ sum(weights);
nulldev = (y' - ybar).^2 * weights;
ka = find(strncmp(gtype,{'covariance','naive'},length(gtype)),1);
if isempty(ka)
    error('unrecognized type');
end

if isempty(offset)
    offset = y * 0;
    is_offset = false;
else
    is_offset = true;
end

if is_sparse
    task = 10;
    [lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(task,parm,x,y-offset,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,ka,cl,intr,maxit,irs,pcs);
else
    task = 11;
    [lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(task,parm,x,y-offset,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,ka,cl,intr,maxit);

end
    
if (jerr ~= 0)
    errmsg = err(jerr,maxit,nx,family);
    if (errmsg.fatal)
        error(errmsg.msg);
    else
        warning(errmsg.msg);
    end
end

ninmax = max(nin);
lam = alm;
if lempty
    lam = fix_lam(lam); % first lambda is infinity; changed to entry point
end

dd=[nvars, lmu];
if ninmax > 0
    ca = ca(1:ninmax,:);
    df = sum(abs(ca) > 0, 1);
    ja = ia(1:ninmax);
    [ja1,oja] = sort(ja);
    beta = zeros(nvars, lmu);
    beta (ja1, :) = ca(oja,:);
else
    beta = zeros(nvars,lmu);
    df = zeros(1,lmu);
end

fit.a0 = a0;
fit.beta = beta;
fit.dev = rsq;
fit.nulldev = nulldev;
fit.df = df';
fit.lambda = lam;
fit.npasses = nlp;
fit.jerr = jerr;
fit.dim = dd;
fit.offset = is_offset;
fit.class = 'elnet';


function new_lam = fix_lam(lam)

new_lam = lam;
if (length(lam) > 2)
    llam=log(lam);
    new_lam(1)=exp(2*llam(2)-llam(3));
end
