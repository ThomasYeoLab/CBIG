function fit = coxnet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
    jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,maxit,family)

% Internal glmnet function. See also glmnet, cvglmnet.

%time --- column 1
%status --- column 2

ty = y(:,1);
tevent = y(:,2);
if (any(ty <= 0))
    error('negative event times encountered;  not permitted for Cox family');
end
if isempty(offset)
    offset = ty * 0;
    is_offset = false;
else
    is_offset = true;
end

if (is_sparse)
    error('Cox model not implemented for sparse x in glmnet');
else
    task = 41;
    [lmu,ca,ia,nin,dev,alm,nlp,jerr,dev0,ot] = glmnetMex(task,parm,x,ty,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,tevent,cl,maxit,offset);

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
    beta (ja1,:) = ca(oja,:);
else
    beta = zeros(nvars,lmu);
    df = zeros(1,lmu);
end

fit.beta = beta;
fit.dev = dev;
fit.nulldev = dev0;
fit.df = df';
fit.lambda = lam;
fit.npasses = nlp;
fit.jerr = jerr;
fit.dim = dd;
fit.offset = is_offset;
fit.class = 'coxnet';


function new_lam = fix_lam(lam)

new_lam = lam;
if (length(lam) > 2)
    llam=log(lam);
    new_lam(1)=exp(2*llam(2)-llam(3));
end

