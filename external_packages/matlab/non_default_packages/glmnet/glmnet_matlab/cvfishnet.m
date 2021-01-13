function result = cvfishnet(object,lambda,x,y,weights,offset,foldid,type,grouped,keep)

% Internal glmnet function. See also cvglmnet.

if nargin < 10 || isempty(keep)
    keep = false;
end

typenames = struct('mse','Mean-Squared Error','mae','Mean Absolute Error','deviance','Poisson Deviance');
if strcmp(type, 'default')
    type = 'deviance';
end
if ~any(strcmp(type, {'mse','mae','deviance'}))
    warning('Only ''mse'', ''deviance'' or ''mae''  available for Poisson models; ''deviance'' used');
    type = 'deviance';
end

is_offset = ~isempty(offset);

predmat = NaN(length(y),length(lambda));
nfolds = max(foldid);
nlams = nfolds;

for i = 1:nfolds
    which = foldid == i;
    fitobj = object{i};
    if (is_offset)
        off_sub = offset(which);
    else
        off_sub = [];
    end
    preds = glmnetPredict(fitobj,x(which,:),[],[],[],off_sub);
    nlami = length(object{i}.lambda);
    predmat(which,1:nlami) = preds;
    nlams(i) = nlami;
end

N = size(y,1) - sum(isnan(predmat),1);

yy = repmat(y, 1, length(lambda));
switch type
    case 'mse'
        cvraw = (yy - predmat).^2;
    case 'mae'
        cvraw = abs(yy - predmat);
    case 'deviance'
        cvraw = devi(yy, predmat);
end

if (length(y)/nfolds < 3) && grouped
    warning('Option grouped=false enforced in cv.glmnet, since < 3 observations per fold');
    grouped = false;
end

if (grouped)
    cvob = cvcompute(cvraw,weights,foldid,nlams);
    cvraw = cvob.cvraw; 
    weights = cvob.weights;
    N = cvob.N;
end

cvm = wtmean(cvraw,weights);
sqccv = (bsxfun(@minus,cvraw,cvm)).^2;
cvsd = sqrt(wtmean(sqccv,weights)./(N-1));

result.cvm = cvm; result.cvsd = cvsd; result.name = typenames.(type);

if (keep)
    result.fit_preval = predmat;
end


function result = devi(yy, eta)

deveta = yy .* eta - exp(eta);
devy = yy .* log(yy) - yy;
devy(yy == 0) = 0;
result = 2 * (devy - deveta);