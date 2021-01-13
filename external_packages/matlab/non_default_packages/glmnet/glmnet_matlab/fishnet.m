function fit = fishnet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
    jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,family)

if any(y < 0)
    error('negative responses encountered;  not permitted for Poisson family');
end

if isempty(offset)
    offset = y * 0;
    is_offset = false;
else
    is_offset = true;
end

if (is_sparse)
    task = 50;
    [lmu,a0,ca,ia,nin,dev,alm,nlp,jerr,dev0,ot] = glmnetMex(task,parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,cl,intr,maxit,offset,irs,pcs);
else
    task = 51;
    [lmu,a0,ca,ia,nin,dev,alm,nlp,jerr,dev0,ot] = glmnetMex(task,parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,cl,intr,maxit,offset);
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
if (ulam == 0.0)
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
fit.dev = dev;
fit.nulldev = dev0;
fit.df = df';
fit.lambda = lam;
fit.npasses = nlp;
fit.jerr = jerr;
fit.dim = dd;
fit.offset = is_offset;
fit.class = 'fishnet';


function new_lam = fix_lam(lam)

new_lam = lam;
if (length(lam) > 2)
    llam=log(lam);
    new_lam(1)=exp(2*llam(2)-llam(3));
end

    
    