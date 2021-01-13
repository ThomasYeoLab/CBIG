function fit = mrelnet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
    jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,jsd,intr,maxit,family)

nr = size(y, 2);

wym = wtmean(y, weights);
nulldev = sum(wtmean(bsxfun(@minus,y,wym).^2,weights)*sum(weights));

if isempty(offset)
    offset = y * 0;
    is_offset = false;
else
    if (size(offset) ~= size(y))
        error('Offset must match dimension of y');
    end
    is_offset = true;
end

if is_sparse
    task = 30;
    [lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(task,parm,x,y-offset,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,cl,intr,maxit,irs,pcs,jsd);
else
    task = 31;
    [lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(task,parm,x,y-offset,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,weights,cl,intr,maxit,jsd);
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

if (nr > 1)
    beta_list = {};
%     a0 = a0 - repmat(mean(a0), nr, 1);  do not center for mrelnet!
    dfmat=a0;
    dd=[nvars, lmu];
    if ninmax > 0
        ca = reshape(ca, nx, nr, lmu);
        ca = ca(1:ninmax,:,:);
        ja = ia(1:ninmax);
        [ja1,oja] = sort(ja);
        df = any(abs(ca) > 0, 2);
        df = sum(df, 1);
        df = df(:);
        for k=1:nr
            ca1 = reshape(ca(:,k,:), ninmax, lmu);
            cak = ca1(oja,:);
            dfmat(k,:) = sum(abs(cak) > 0, 1);
            beta = zeros(nvars, lmu);
            beta(ja1,:) = cak;
            beta_list{k} = beta;
        end
    else
        for k = 1:nr
            dfmat(k,:) = zeros(1,lmu);
            beta_list{k} = zeros(nvars, lmu);
        end
        df = zeros(1,lmu);
    end
    fit.beta = beta_list;
    fit.dfmat = dfmat;
    
    
else
    
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
    
    fit.beta = beta;
    
    
end

fit.a0 = a0;
fit.dev = rsq;
fit.nulldev = nulldev;
fit.df = df';
fit.lambda = lam;
fit.npasses = nlp;
fit.jerr = jerr;
fit.dim = dd;
fit.offset = is_offset;
fit.class = 'mrelnet';


function new_lam = fix_lam(lam)

new_lam = lam;
if (length(lam) > 2)
    llam=log(lam);
    new_lam(1)=exp(2*llam(2)-llam(3));
end

