function result = cvelnet(object, lambda, x, y, weights, offset, foldid, ...
    type, grouped, keep)

% Internal glmnet function. See also cvglmnet.

if nargin < 10 || isempty(keep)
    keep = false;
end

typenames = struct('deviance','Mean-Squared Error','mse','Mean-Squared Error','mae','Mean Absolute Error');
if strcmp(type, 'default')
    type = 'mse';
end
if ~any(strcmp(type, {'mse','mae','deviance'}))
    warning('Only ''mse'', ''deviance'' or ''mae''  available for Gaussian models; ''mse'' used');
    type = 'mse';
end

if ~isempty(offset)
    y = y - offset;
end

predmat = NaN(length(y),length(lambda));
nfolds = max(foldid);
nlams = nfolds;

for i = 1:nfolds
    which = foldid == i;
    fitobj = object{i};
    fitobj.offset = false;
    preds = glmnetPredict(fitobj,x(which,:));
    nlami = length(object{i}.lambda);
    predmat(which,1:nlami) = preds;
    nlams(i) = nlami;
end

N = size(y,1) - sum(isnan(predmat),1);

yy = repmat(y, 1, length(lambda));
switch type
    case 'mse'
        cvraw = (yy - predmat).^2;
    case 'deviance'
        cvraw = (yy - predmat).^2;
    case 'mae'
        cvraw = abs(yy - predmat);
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
end