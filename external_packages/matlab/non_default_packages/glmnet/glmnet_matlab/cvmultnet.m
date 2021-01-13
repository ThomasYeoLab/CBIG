function result = cvmultnet(object, lambda, x, y, weights, offset, foldid, ...
    type, grouped, keep)

if nargin < 10 || isempty(keep)
    keep = false;
end

typenames = struct('mse','Mean-Squared Error','mae','Mean Absolute Error',...
    'deviance','Multinomial Deviance','class','Misclassification Error');

if strcmp(type,'default')
    type = 'deviance';
end

if ~any(strcmp(type,{'mse','mae','deviance','class'}))
    warning('Only ''deviance'', ''class'', ''mse'' or ''mae''  available for multinomial models; ''deviance'' used');
    type = 'deviance';
end

prob_min = 1e-5; prob_max = 1 - prob_min;
nc = size(y);
if nc(2) == 1
    [classes,~,sy] = unique(y);
    nc = length(classes);
    indexes = eye(nc);
    y = indexes(sy,:);
else
    nc = nc(2);
end
is_offset = ~isempty(offset);
predmat = NaN(size(y,1),nc,length(lambda));
nfolds = max(foldid);
nlams = zeros(nfolds,1);

for i = 1:nfolds
    which = foldid==i;
    fitobj = object{i};
    if (is_offset)
        off_sub = offset(which,:);
    else
        off_sub = [];
    end
    preds = glmnetPredict(fitobj,x(which,:),[],'response',[],off_sub);
    nlami = length(object{i}.lambda);
    predmat(which,:,1:nlami) = preds;
    nlams(i) = nlami;
end

ywt = sum(y, 2);
y = y ./ repmat(ywt,1,size(y,2));
weights = weights .* ywt;
N = size(y,1) - sum(isnan(predmat(:,1,:)),1);
bigY = repmat(y, [1,1,length(lambda)]);
switch type
    case 'mse'
        cvraw = squeeze(sum((bigY - predmat).^2, 2));
    case 'mae'
        cvraw = squeeze(sum(abs(bigY - predmat), 2));
    case 'deviance'
        predmat = min(max(predmat,prob_min),prob_max);
        lp = bigY .* log(predmat);
        ly = bigY .* log(bigY);
        ly(bigY == 0) = 0;
        cvraw = squeeze(sum(2 * (ly - lp), 2));
    case 'class'
        classid = NaN(size(y,1),length(lambda));
        for i = 1:length(lambda)
            classid(:,i) = glmnet_softmax(predmat(:,:,i));
        end
        classid = reshape(classid,[],1);
        yperm = reshape(permute(bigY, [1,3,2]),[],nc);
        idx = sub2ind(size(yperm), 1:length(classid), classid');
        cvraw = reshape(1 - yperm(idx), [], length(lambda));
end

if (grouped)
    cvob = cvcompute(cvraw, weights, foldid, nlams);
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
        
        



    


    