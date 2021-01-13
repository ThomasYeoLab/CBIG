function fit = lognet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
    jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,kopt,family)

[noo,nc] = size(y);
if noo ~= nobs
    error('x and y have different number of rows in call to glmnet');
end
if nc == 1
    [classes,~,sy] = unique(y);
    nc = length(classes);
    indexes = eye(nc);
    y = indexes(sy,:);
else
    classes = 1: nc;
end

if strcmp(family, 'binomial')
    if nc > 2
        error ('More than two classes; use multinomial family instead');
    end
    nc = 1; % for calling binomial
    y = y(:,[2,1]);
end
o = [];
if ~isempty(weights)
    % check if any are zero
    o = weights > 0;
    if ~all(o) %subset the data
        y = y(o,:);
        x = x(o,:);
        weights = weights(o);
        nobs = sum(o);
    else
        o = [];
    end
    [my,ny] = size(y);
    y = y .* repmat(weights,1,ny);
end

if isempty(offset)
    offset = y * 0;
    is_offset = false;
else
    if ~isempty(o)
        offset = offset(o,:);
    end
    do = size(offset);
    if (do(1) ~= nobs)
        error('offset should have the same number of values as observations in binomial/multinomial call to glmnet');
    end
    if (nc == 1)
        if (do(2) == 1)
            offset = cat(2,offset,-offset);
        end
        if (do(2) > 2)
            error('offset should have 1 or 2 columns in binomial call to glmnet');
        end
    end
    if strcmp(family,'multinomial') && (do(2) ~= nc)
        error('offset should have same shape as y in multinomial call to glmnet');
    end
    is_offset = true;
end

if (is_sparse)
    task = 20;
    [lmu,a0,ca,ia,nin,dev,alm,nlp,jerr,dev0,ot] = glmnetMex(task,parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,cl,intr,maxit,nc,kopt,offset,irs,pcs);
else
    task = 21;
    [lmu,a0,ca,ia,nin,dev,alm,nlp,jerr,dev0,ot] = glmnetMex(task,parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thresh,isd,cl,intr,maxit,nc,kopt,offset);
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

if strcmp(family, 'multinomial')
    beta_list = {};
    a0 = a0 - repmat(mean(a0), nc, 1);      %multinomial: center the coefficients
    dfmat=a0;
    dd=[nvars, lmu];
    if ninmax > 0
        ca = reshape(ca, nx, nc, lmu);
        ca = ca(1:ninmax,:,:);
        ja = ia(1:ninmax);
        [ja1,oja] = sort(ja);
        df = any(abs(ca) > 0, 2);
        df = sum(df, 1);
        df = df(:)';
        for k=1:nc
            ca1 = reshape(ca(:,k,:), ninmax, lmu);
            cak = ca1(oja,:);
            dfmat(k,:) = sum(abs(cak) > 0, 1);
            beta = zeros(nvars, lmu);
            beta(ja1,:) = cak;
            beta_list{k} = beta;
        end
    else
        for k = 1:nc
            dfmat(k,:) = zeros(1,lmu);
            beta_list{k} = zeros(nvars, lmu);
        end
        df = zeros(1,lmu);
    end
    fit.a0 = a0;
    fit.label = classes;
    fit.beta = beta_list;
    fit.dev = dev;
    fit.nulldev = dev0;
    fit.dfmat = dfmat;
    fit.df = df';
    fit.lambda = lam;
    fit.npasses = nlp;
    fit.jerr = jerr;
    fit.dim = dd;
    if (kopt == 2)
        grouped = true;
    else
        grouped = false;
    end
    fit.grouped = grouped;
    fit.offset = is_offset;
    fit.class = 'multnet';
    
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
    fit.a0 = a0;
    fit.label = classes;
    fit.beta = beta; %sign flips make 2 arget class
    fit.dev = dev;
    fit.nulldev = dev0;
    fit.df = df';
    fit.lambda = lam;
    fit.npasses = nlp;
    fit.jerr = jerr;
    fit.dim = dd;
    fit.offset = is_offset;
    fit.class = 'lognet';
end


function new_lam = fix_lam(lam)

new_lam = lam;
if (length(lam) > 2)
    llam=log(lam);
    new_lam(1)=exp(2*llam(2)-llam(3));
end
