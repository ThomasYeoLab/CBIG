function result = cvmrelnet(object, lambda, x, y, weights, offset, foldid, ...
    type, grouped, keep)

if nargin < 10 || isempty(keep)
    keep = false;
end

typenames = struct('deviance','Mean-Squared Error','mse','Mean-Squared Error',...
    'mae','Mean Absolute Error');

if strcmp(type,'default')
    type = 'mse';
end

if ~any(strcmp(type,{'mse','mae','deviance'}))
    warning('Only ''mse'',''deviance'' or ''mae'' available for multiresponse Gaussian models; ''mse'' used');
    type = 'mse';
end

[nobs, nc] = size(y);

if ~isempty(offset)
    y = y - offset;
end

predmat = NaN(nobs,nc,length(lambda));
nfolds = max(foldid);
nlams = nfolds;

for i = 1:nfolds
    which = foldid == i;
    fitobj = object{i};
    fitobj.offset = false;
    preds = glmnetPredict(fitobj,x(which,:));
    nlami = length(object{i}.lambda);
    predmat(which,:,1:nlami) = preds;
    nlams(i) = nlami;
end
N = nobs - reshape(sum(isnan(predmat(:,1,:)),1),1,[]);
bigY = repmat(y, [1,1,length(lambda)]);
switch type
    case 'mse'
        cvraw = squeeze(sum((bigY - predmat).^2, 2));
    case 'mae'
        cvraw = squeeze(sum(abs(bigY - predmat), 2));
end
if (nobs/nfolds < 3) && grouped
    warning('Option grouped=false enforced in cv.glmnet, since < 3 observations per fold');
    grouped = false;
end
if (grouped)
 cvob = cvcompute(cvraw, weights, foldid, nlams);
    cvraw = cvob.cvraw;
    weights = cvob.weights;
    N = cvob.N;
end
% end
cvm = wtmean(cvraw,weights);
sqccv = (bsxfun(@minus,cvraw,cvm)).^2;
cvsd = sqrt(wtmean(sqccv,weights)./(N-1));
result.cvm = cvm; result.cvsd = cvsd; result.name = typenames.(type);
if (keep)
    result.fit_preval = predmat;
end   


function result = glmnet_softmax(x)

d = size(x);
nas = any(isnan(x),2);
if any(nas)
    pclass = NaN(d(1),1);
    if (sum(nas) < d(1))
        pclass2 = glmnet_softmax(x(~nas,:));
        pclass(~nas) = pclass2;
        result = pclass;
    end
else
    maxdist = x(:,1);
    pclass = ones(d(1),1);
    for i = 2:d(2)
        l = x(:,i)>maxdist;
        pclass(l) = i;
        maxdist(l) = x(l,i);
    end
    result = pclass;
end



