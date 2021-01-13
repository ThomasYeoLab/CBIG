function result = cvcoxnet(object, lambda, x, y, weights, offset, foldid, ...
    type, grouped, keep)

% Internal glmnet function. See also cvglmnet.

if nargin < 10 || isempty(keep)
    keep = false;
end

typenames = struct('deviance','Partial Likelihood Deviance');
if strcmp(type, 'default')
    type = 'deviance';
end
if ~any(strcmp(type, {'deviance'}))
    warning('Only ''deviance''  available for Cox models; changed to type=''deviance''');
    type = 'deviance';
end

if isempty(offset)
    offset = zeros(size(x,1),1);
end
nfolds = max(foldid);

if (length(weights)/nfolds < 10) && ~grouped
    warning('Option grouped=true enforced for cv.coxnet, since < 3 observations per fold');
    grouped = true;
end
cvraw = NaN(nfolds,length(lambda));

for i = 1:nfolds
    which = foldid == i;
    fitobj = object{i};
    coefmat = glmnetPredict(fitobj,[],[],'coefficients');
    if (grouped)
        plfull = cox_deviance([],y,x,offset,weights,coefmat);
        plminusk = cox_deviance([],y(~which,:),x(~which,:),offset(~which),...
            weights(~which),coefmat);
        cvraw(i,1:length(plfull)) = plfull - plminusk;
    else
        plk = cox_deviance([],y(which,:),x(which,:),offset(which),...
            weights(which),coefmat);
        cvraw(i,1:length(plk)) = plk;
    end
end
status = y(:,2);
N = nfolds - sum(isnan(cvraw),1);
weights = accumarray(reshape(foldid,[],1),weights.*status);
cvraw = bsxfun(@rdivide,cvraw,weights);  %even some weight = 0 does matter because of adjustment in wtmean!
cvm = wtmean(cvraw,weights);
sqccv = (bsxfun(@minus,cvraw,cvm)).^2;
cvsd = sqrt(wtmean(sqccv,weights)./(N-1));
result.cvm = cvm; result.cvsd = cvsd; result.name = typenames.(type);
if (keep)
    warning('keep=TRUE not implemented for coxnet');
end


function result = cox_deviance(pred, y, x, offset, weights, beta)

ty = y(:,1);
tevent = y(:,2);
nobs = length(ty); nvars = size(x,2);
if isempty(weights)
    weights = ones(nobs,1);
end
if isempty(offset)
    offset = zeros(nobs,1);
end
if isempty(beta)
    beta = []; nvec = 1; nvars = 0;
else
    nvec = size(beta,2);
end

task = 42;
[flog, jerr] = glmnetMex(task,x,ty,tevent,offset,weights,nvec,beta);


if (jerr ~= 0)
    errmsg = err(jerr,0,0,'cox');
    if (errmsg.fatal)
        error(errmsg.msg);
    else
        warning(errmsg.msg);
    end
end
result = -2 * flog;