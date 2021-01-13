function CVerr = cvglmnet(x,y,family,options,type,nfolds,foldid,parallel,keep,grouped)

%--------------------------------------------------------------------------
% cvglmnet.m: cross-validation for glmnet
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Does k-fold cross-validation for glmnet, produces a plot, and returns
%    a value for lambda.
%
% USAGE:
%    CVerr = cvglmnet(x, y, family, options, type, nfolds, foldid,
%    parallel, keep, grouped);
%
%    Fewer input arguments(more often) are allowed in the call, but must
%    come in the order listed above. To set default values on the way, use
%    empty matrix []. 
%    For example, CVfit=cvglmnet(x,y,'multinomial',[],[],20).
%    
%
% INPUT ARGUMENTS
% x           x matrix as in glmnet.
% y           Response y as in glmnet.
% family      Response type as family in glmnet.
% options     Options as in glmnet.
% type        loss to use for cross-validation. Currently five options, not
%             all available for all models. The default is type='deviance', which uses
%             squared-error for Gaussian models (a.k.a type='mse' there), deviance for
%             logistic and Poisson regression, and partial-likelihood for the Cox
%             model. type='class' applies to binomial and multinomial logistic
%             regression only, and gives misclassification error. type='auc' is for
%             two-class logistic regression only, and gives area under the ROC curve.
%             type='mse' or type='mae' (mean absolute error) can be used by all models
%             except the 'cox'; they measure the deviation from the fitted mean to the
%             response.  
% nfolds      number of folds - default is 10. Although nfolds can be as
%             large as the sample size (leave-one-out CV), it is not recommended for
%             large datasets. Smallest value allowable is nfolds=3.
% foldid      an optional vector of values between 1 and nfold identifying
%             what fold each observation is in. If supplied, nfold can be
%             missing.
% parallel    If true, use parallel computation to fit each fold. If a
%             worker pool is not open, it will open using the
%             default cluster profile and close after the computation is
%             over. 
% keep        If keep=true, a prevalidated array is returned containing
%             fitted values for each observation and each value of lambda.
%             This means these fits are computed with this observation and
%             the rest of its fold omitted. The foldid vector is also
%             returned. Default is keep=false.   
% grouped     This is an experimental argument, with default true, and can
%             be ignored by most users. For all models except the 'cox',
%             this refers to computing nfolds separate statistics, and then
%             using their mean and estimated standard error to describe the
%             CV curve. If grouped=false, an error matrix is built up at
%             the observation level from the predictions from the nfold
%             fits, and then summarized (does not apply to
%             type='auc'). For the 'cox' family, grouped=true obtains the 
%             CV partial likelihood for the Kth fold by subtraction; by
%             subtracting the log partial likelihood evaluated on the full
%             dataset from that evaluated on the on the (K-1)/K dataset.
%             This makes more efficient use of risk sets. With
%             grouped=FALSE the log partial likelihood is computed only on
%             the Kth fold.
%
% OUTPUT ARGUMENTS:
% A structure is returned with the following fields.
% lambda      the values of lambda used in the fits.
% cvm         the mean cross-validated error - a vector of length
%             length(lambda). 
% cvsd        estimate of standard error of cvm.
% cvup        upper curve = cvm+cvsd.
% cvlo        lower curve = cvm-cvsd.
% nzero       number of non-zero coefficients at each lambda.
% name        a text string indicating type of measure (for plotting
%             purposes). 
% glmnet_fit  a fitted glmnet object for the full data.
% lambda_min  value of lambda that gives minimum cvm.
% lambda_1se  largest value of lambda such that error is within 1 standard
%             error of the minimum. 
% class       Type of regression - internal usage.
% fit_preval  if keep=true, this is the array of prevalidated fits. Some
%             entries can be NA, if that and subsequent values of lambda
%             are not reached for that fold.
% foldid      if keep=true, the fold assignments used.
%
% DETAILS:
%    The function runs glmnet nfolds+1 times; the first to get the lambda
%    sequence, and then the remainder to compute the fit with each of the 
%    folds omitted. The error is accumulated, and the average error and 
%    standard deviation over the folds is computed. Note that cv.glmnet 
%    does NOT search for values for alpha. A specific value should be 
%    supplied, else alpha=1 is assumed by default. If users would like to 
%    cross-validate alpha as well, they should call cv.glmnet with a 
%    pre-computed vector foldid, and then use this same fold vector in 
%    separate calls to cv.glmnet with different values of alpha. 
%
% LICENSE: GPL-2
%
% DATE: 30 Aug 2013
%
% AUTHORS:
%    Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani
%    Fortran code was written by Jerome Friedman
%    R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
%    The original MATLAB wrapper was written by Hui Jiang (14 Jul 2009),
%    and was updated and is maintained by Junyang Qian (30 Aug 2013) junyangq@stanford.edu,
%    Department of Statistics, Stanford University, Stanford, California, USA.
%
% REFERENCES:
%    Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, 
%    http://www.jstatsoft.org/v33/i01/
%    Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
%    
%    Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent,
%    http://www.jstatsoft.org/v39/i05/
%    Journal of Statistical Software, Vol. 39(5) 1-13
%
%    Tibshirani, Robert., Bien, J., Friedman, J.,Hastie, T.,Simon, N.,Taylor, J. and Tibshirani, Ryan. (2010) Strong Rules for Discarding Predictors in Lasso-type Problems,
%    http://www-stat.stanford.edu/~tibs/ftp/strong.pdf
%    Stanford Statistics Technical Report
%
% SEE ALSO:
%    cvglmnetPlot, cvglmnetCoef, cvglmnetPredict, and glmnet.
%
% EXAMPLES:
%    n=1000; p=100;
%    nzc=fix(p/10);
%    x=randn(n,p);
%    beta=randn(nzc,1);
%    fx=x(:,1:nzc) * beta;
%    eps=randn(n,1)*5;
%    y=fx+eps;
%    px=exp(fx);
%    px=px./(1+px);
%    ly=binornd(1,px,length(px),1);   
%    cvob1=cvglmnet(x,y);
%    cvglmnetPlot(cvob1);
%    cvglmnetCoef(cvob1)
%    cvglmnetPredict(cvob1,x(1:5,:),'lambda_min')
% 
%    cvobla=cvglmnet(x,y,[],[],'mae');
%    cvglmnetPlot(cvobla);
%    
%    cvob2=cvglmnet(x,ly,'binomial');
%    cvglmnetPlot(cvob2);
%    
%    figure;
%    cvob3=cvglmnet(x,ly,'binomial',[],'class');
%    cvglmnetPlot(cvob3);
% 
%    mu=exp(fx/10);
%    y=poissrnd(mu,n,1);
%    cvob4=cvglmnet(x,y,'poisson');
%    cvglmnetPlot(cvob4);
%    
% % Multinomial
%    n=500; p=30;
%    nzc=fix(p/10);
%    x=randn(n,p);
%    beta3=randn(10,3);
%    beta3=cat(1,beta3,zeros(p-10,3));
%    f3=x*beta3;
%    p3=exp(f3);
%    p3=bsxfun(@rdivide,p3,sum(p3,2));
%    g3=mnrnd(1,p3);
%    g3=g3*(1:size(p3,2))';
%    cvfit=cvglmnet(x,g3,'multinomial');
%    cvglmnetPlot(cvfit);
%    
% % Cox
%    n=1000;p=30;
%    nzc=p/3;
%    x=randn(n,p);
%    beta=randn(nzc,1);
%    fx=x(:,1:nzc)*beta/3;
%    hx=exp(fx);
%    ty=exprnd(1./hx,n,1);
%    tcens=binornd(1,0.3,n,1);
%    y=cat(2,ty,1-tcens);
%    foldid=randsample(10,n,true);
%    fit1_cv=cvglmnet(x,y,'cox',[],[],[],foldid);
%    cvglmnetPlot(fit1_cv);
%    
% % Parallel
%    matlabpool;
%    x=randn(1e3,100);
%    y=randn(1e3,1);
%    tic;
%    cvglmnet(x,y);
%    toc;
%    tic;
%    cvglmnet(x,y,[],[],[],[],[],true);
%    toc;
%
% DEVELOPMENT:
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported.
%    29 Dec 2013: Fixed a bug in the return value of CVerr.fit_preval,
%                 pointed out by Leon Peshkin from Harvard University.
%
% OLDER UPDATES:
%    26 Jan 2010: Fixed a bug in the description of y, pointed out by
%                 Peter Rijnbeek from Erasmus University.
%    09 Mar 2010: Fixed a bug of printing "ka = 2", pointed out by
%                 Ramon Casanova from Wake Forest University.
%    25 Mar 2010: Fixed a bug when p > n in multinomial fitting, pointed
%                 out by Gerald Quon from University of Toronto
%    25 Jul 2010: Check for input matrix format and size
%    27 Sep 2010: Fixed a bug of undefined "df" in multinomial fitting,
%                 pointed by Jeff Howbert from Insilicos.

%%% Set default values
if nargin < 3 || isempty(family)
    family = 'gaussian';
end
if nargin < 4
    options = [];    
end
if nargin < 5 || isempty(type)
    type = 'default';
end
if nargin < 6 || isempty(nfolds)
    nfolds = 10;
end
if nargin < 7
    foldid = [];
end
if nargin < 8 || isempty(parallel)
    parallel = false;
end
if nargin < 9 || isempty(keep)
    keep = false;
end
if nargin < 10 || isempty(grouped)
    grouped = true;
end

options = glmnetSet(options);

if (~isempty(options.lambda)) && (length(options.lambda)<2)
    error('Need more than one value of lambda for cv.glmnet');
end

N = size(x, 1);

if (size(y,1) ~= N)
    y = transpose(y);
end

if (~isempty(options.offset)) && (size(options.offset, 1) ~= N)
    options.offset = transpose(options.offset);
end

if (isempty(options.weights))
    options.weights = ones(N,1);
end

glmfit = glmnet(x, y, family, options);

is_offset = glmfit.offset;
options.lambda = glmfit.lambda;

nz = glmnetPredict(glmfit,[],[],'nonzero');
if (strcmp(glmfit.class,'multnet'))
    nnz = zeros(length(options.lambda),length(nz));
    for i = 1:length(nz)
        nnz(:,i) = transpose(sum(nz{i},1));
    end
    nz = ceil(median(nnz,2));
elseif strcmp(glmfit.class, 'mrelnet') 
    nz = transpose(sum(nz{1}, 1));
else
    nz = transpose(sum(nz,1));
end

if isempty(foldid)
    population = cat(2, repmat(1:nfolds, 1, floor(N/nfolds)), 1:mod(N,nfolds));
    foldid = population(randperm(length(population), N));
else
    nfolds = max(foldid);
end
foldid = reshape(foldid, numel(foldid), 1);

if (nfolds < 3)
    error('nfolds must be bigger than 3; nfolds=10 recommended');
end

cpredmat = cell(nfolds,1);

if (parallel == true)
    offpar = 0;
    if matlabpool('size') <= 0
        offpar = 1;
        matlabpool;
    end
    
    parfor i = 1: nfolds
        which = foldid==i;
        opts = options;
        opts.weights = opts.weights(~which,:);
        opts.lambda = options.lambda;
        if (is_offset)
            opts.offset = opts.offset(~which,:);
        end
        xr = x(~which,:); yr = y(~which,:);
        cpredmat{i} = glmnet(xr, yr, family, opts);
    end
    
    if (offpar)
        matlabpool close;
    end    
else   
    for i = 1: nfolds        
        which = foldid==i;
        opts = options;
        opts.weights = opts.weights(~which,:);
        opts.lambda = options.lambda;
        if (is_offset)
            opts.offset = opts.offset(~which,:);
        end
        xr = x(~which,:); yr = y(~which,:);
        cpredmat{i} = glmnet(xr, yr, family, opts);
    end
end

switch cpredmat{1}.class
    case 'elnet'
        cvstuff = cvelnet(cpredmat,options.lambda,x,y,options.weights,options.offset,foldid,type,grouped,keep);
    case 'lognet'
        cvstuff = cvlognet(cpredmat,options.lambda,x,y,options.weights,options.offset,foldid,type,grouped,keep);
    case 'multnet'
        cvstuff = cvmultnet(cpredmat,options.lambda,x,y,options.weights,options.offset,foldid,type,grouped,keep);
    case 'coxnet'
        cvstuff = cvcoxnet(cpredmat,options.lambda,x,y,options.weights,options.offset,foldid,type,grouped,keep);
    case 'mrelnet'
        cvstuff = cvmrelnet(cpredmat,options.lambda,x,y,options.weights,options.offset,foldid,type,grouped,keep); 
    case 'fishnet'
        cvstuff = cvfishnet(cpredmat,options.lambda,x,y,options.weights,options.offset,foldid,type,grouped,keep);
end

cvm = cvstuff.cvm;
cvsd = cvstuff.cvsd;
cvname = cvstuff.name;

CVerr.lambda = options.lambda;
CVerr.cvm = transpose(cvm); CVerr.cvsd = transpose(cvsd); 
CVerr.cvup = transpose(cvm+cvsd); CVerr.cvlo = transpose(cvm-cvsd); CVerr.nzero = nz;
CVerr.name = cvname; CVerr.glmnet_fit = glmfit;
if (keep)
    CVerr.fit_preval = cvstuff.fit_preval; CVerr.foldid = foldid;
end
if strcmp(type, 'auc')
    cvm = -cvm;
end
CVerr.lambda_min = max(options.lambda(cvm<=min(cvm)));
idmin = options.lambda==CVerr.lambda_min;
semin = cvm(idmin)+cvsd(idmin);
CVerr.lambda_1se = max(options.lambda(cvm<=semin));
CVerr.class = 'cv.glmnet';
end



