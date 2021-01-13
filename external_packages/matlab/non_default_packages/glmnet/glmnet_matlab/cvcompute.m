function cvcpt = cvcompute(mat, weights, foldid, nlams)

% Internal glmnet function. See also cvglmnet.

% Compute the weighted mean and SD within folds, and hence the se of the mean

wisum = accumarray(foldid,weights);
nfolds = max(foldid);
outmat = NaN(nfolds,size(mat,2));
good = zeros(nfolds,size(mat,2));
mat(isinf(mat)) = NaN;
for i = 1:nfolds
    mati = mat(foldid==i,:);
    wi = weights(foldid==i,:);
    outmat(i,:) = wtmean(mati,wi);
    good(i,1:nlams(i)) = 1;
end
N = sum(good,1);
cvcpt.cvraw = outmat;
cvcpt.weights = wisum;
cvcpt.N = N;
end