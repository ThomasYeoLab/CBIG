function fit = glmnet(x, y, family, options)

%--------------------------------------------------------------------------
% glmnet.m: fit an GLM with lasso or elasticnet regularization
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Fit a generalized linear model via penalized maximum likelihood. The 
%    regularization path is computed for the lasso or elasticnet penalty 
%    at a grid of values for the regularization parameter lambda. Can deal 
%    with all shapes of data, including very large sparse data matrices. 
%    Fits linear, logistic and multinomial, Poisson, and Cox regression 
%    models.
%
% USAGE:
%    fit = glmnet(x, y)
%    fit = glmnet(x, y, family)
%    fit = glmnet(x, y, family, options)
%    (Use empty matrix [] to apply the default value, eg. fit = glmnet(x,
%    y, [], options))
%
%
% EXTERNAL FUNCTIONS:
% options         = glmnetSet;                  provided with glmnet.m
%
% INPUT ARGUMENTS:
% x           Input matrix, of dimension nobs x nvars; each row is an
%             observation vector. Can be in sparse matrix format.
% y           Response variable. Quantitative (column vector) for family =
%             'gaussian', or family = 'poisson'(non-negative counts). For
%             family = 'binomial' should be either a column vector with two 
%             levels, or a two-column matrix of counts or proportions. For
%             family = 'multinomial', can be a column vector of nc>=2
%             levels, or a matrix with nc columns of counts or proportions.
%             For family = 'cox', y should be a two-column matrix with the
%             first column for time and the second for status. The latter
%             is a binary variable, with 1 indicating death, and 0
%             indicating right censored. For family = 'mgaussian', y is a
%             matrix of quantitative responses. 
% family      Reponse type. (See above). Default is 'gaussian'.
% options     A structure that may be set and altered by glmnetSet.
%             Default values for some often used options:
%                options.alpha = 1.0  (elastic-net mixing parameter)
%                options.nlambda = 100  (number of lambda values)
%                options.lambda depends on data, nlambda and
%                lambda_min(user spplied lambda sequence) 
%                options.standardize = true  (variable standardization)
%                options.weights = all ones vector (observation weights)
%             For more details, type help glmnetSet.
%             
% OUTPUT ARGUMENTS:
% fit         A structure.
% fit.a0      Intercept sequence of length length(fit.lambda).
% fit.beta    For "elnet" and "lognet" models, a nvars x length(lambda)
%             matrix of coefficients. For "multnet", a list of nc such
%             matrices, one for each class.
% fit.lambda  The actual sequence of lambda values used.
% fit.dev     The fraction of (null) deviance explained (for "elnet", this
%             is the R-square).
% fit.nulldev Null deviance (per observation).
% fit.df      The number of nonzero coefficients for each value of lambda.
%             For "multnet", this is the number of variables with a nonzero
%             coefficient for any class.
% fit.dfmat   For "multnet" only. A matrix consisting of the number of
%             nonzero coefficients per class.
% fit.dim     Dimension of coefficient matrix (ices).
% fit.npasses Total passes over the data summed over all lambda values.
% fit.offset  a logical variable indicating whether an offset was included
%             in the model.
% fit.jerr    Error flag, for warnings and errors (largely for internal
%             debugging).
% fit.class   Type of regression - internal usage.
% fit.call    a cell including the names of all the input variables in the
%             parent environment.
%
% DETAILS:
%    The sequence of models implied by lambda is fit by coordinate descent.
%    For family='gaussian' this is the lasso sequence if alpha=1, else it
%    is the elasticnet sequence. For the other families, this is a lasso or
%    elasticnet regularization path for fitting the generalized linear
%    regression paths, by maximizing the appropriate penalized
%    log-likelihood (partial likelihood for the 'cox' model). Sometimes the
%    sequence is truncated before nlambda values of lambda have been used,
%    because of instabilities in the inverse link functions near a
%    saturated fit. glmnet(...,family='binomial') fits a traditional
%    logistic regression model for the log-odds.
%    glmnet(...,family='multinomial') fits a symmetric multinomial model,
%    where each class is represented by a linear model (on the log-scale).
%    The penalties take care of redundancies. A two-class 'multinomial'
%    model will produce the same fit as the corresponding 'binomial' model,
%    except the pair of coefficient matrices will be equal in magnitude and
%    opposite in sign, and half the 'binomial' values. Note that the
%    objective function for 'gaussian' is
%
%                    1/2 RSS / nobs + lambda * penalty,
%
%    and for the logistic models it is
%
%                    -loglik / nobs + lambda * penalty.
%
%    Note also that for 'gaussian', glmnet standardizes y to have unit
%    variance before computing its lambda sequence (and then unstandardizes
%    the resulting coefficients); if you wish to reproduce/compare results
%    with other software, best to supply a standardized y. The latest two
%    features in glmnet are the family='mgaussian' family and the
%    mtype='grouped' in options for multinomial fitting. The former
%    allows a multi-response gaussian model to be fit, using a "group
%    -lasso" penalty on the coefficients for each variable. Tying the
%    responses together like this is called "multi-task" learning in some
%    domains. The grouped multinomial allows the same penalty for the
%    family='multinomial' model, which is also multi-responsed. For both of
%    these the penalty on the coefficient vector for variable j is
%
%            (1-alpha)/2 * ||beta_j||_2^2 + alpha * ||beta_j||_2
%
%    When alpha=1 this is a group-lasso penalty, and otherwise it mixes
%    with quadratic just like elasticnet. 
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
%    glmnetPrint, glmnetPlot, glmnetCoef, glmnetPredict,
%    glmnetSet, glmnetControl and cvglmnet.
%
% EXAMPLES:
% % Gaussian
%    x=randn(100,20);
%    y=randn(100,1);
%    fit1 = glmnet(x,y);
%    glmnetPrint(fit1);
%    glmnetPredict(fit1,[],0.01,'coef')  %extract coefficients at a single value of lambda
%    glmnetPredict(fit1,x(1:10,:),[0.01,0.005]')  %make predictions
%
% % Multivariate Gaussian:
%    y=randn(100,3);
%    fit1m=glmnet(x,y,'mgaussian');
%    glmnetPlot(fit1m,[],[],'2norm');
%
% % Binomial:
%    g2=randsample(2,100,true);
%    fit2=glmnet(x,g2,'binomial');
%
% % Multinomial:
%    g4=randsample(4,100,true);
%    fit3=glmnet(x,g4,'multinomial');
%    opts=struct('mtype','grouped');
%    fit3a=glmnet(x,g4,'multinomial',opts);
%
% % Poisson:
%    N=500; p=20;
%    nzc=5;
%    x=randn(N,p);
%    beta=randn(nzc,1);
%    f=x(:,1:nzc) * beta;
%    mu=exp(f);
%    y=poissrnd(mu,N,1);
%    fit=glmnet(x,y,'poisson');
%    glmnetPlot(fit);
%    pfit=glmnetPredict(fit,x,0.001,'response');
%    plot(pfit,y,'o');
%
% % Cox:
%    N=1000; p=30;
%    nzc=p/3;
%    x=randn(N,p);
%    beta=randn(nzc,1);
%    fx=x(:,1:nzc)*beta/3;
%    hx=exp(fx);
%    ty=exprnd(1./hx,N,1);
%    tcens=binornd(1,0.3,N,1);
%    y=cat(2,ty,1-tcens);
%    fit=glmnet(x,y,'cox');
%    glmnetPlot(fit);
%
% % Sparse:
%    n=10000;p=200;
%    nzc=fix(p/10);
%    x=randn(n,p);
%    iz=randsample(n*p,n*p*0.85,false);
%    x(iz)=0;
%    sx=sparse(x);
%    beta=randn(nzc,1);
%    fx=x(:,1:nzc)*beta;
%    eps=randn(n,1);
%    y=fx+eps;
%    px=exp(fx);
%    px=px./(1+px);
%    ly=binornd(1,px,length(px),1);
%    tic;
%    fit1=glmnet(sx,y);
%    toc;
%    tic;
%    fit2n=glmnet(x,y);
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

if nargin < 2
    error('more input arguments needed.');
end

if nargin < 3 || isempty(family)
    family = 'gaussian';
end

if nargin < 4 || isempty(options)
    options = glmnetSet;
end

%Get the names of input variables
out_x = inputname(1); out_y = inputname(2);
out_family = mat2str([]); out_options = mat2str([]);
if nargin > 2
    if ~isempty(inputname(3))
        out_family = inputname(3);
    else
        out_family = family;
    end
end
if nargin > 3
    if ~isempty(inputname(4))
        out_options = inputname(4);
    end
end

%match the family, abbreviation allowed
fambase = {'gaussian','binomial','poisson','multinomial','cox','mgaussian'};
famind = find(strncmp(family,fambase,length(family)),1);
if isempty(famind)
    error('family should be one of ''gaussian'', ''binomial'', ''poisson'', ''multinomial'', ''cox'', ''mgaussian''');
else
    family = fambase{famind};
end

% Prepare parameters
options = glmnetSet(options);

if (options.alpha > 1)
    warning('alpha >1; set to 1');
    options.alpha = 1;
end
if (options.alpha < 0)
    warning('alpha <0; set to 0');
    options.alpha = 0;
end

parm = options.alpha;
nlam = options.nlambda;
[nobs,nvars] = size(x);

weights = options.weights;
if isempty(weights)
    weights = ones(nobs,1);
else
    if (length(weights) ~= nobs)
        error('number of elements in weights (%d) not equal to the number of rows of x (%d)',length(weights),nobs);
    end
end

nrowy = size(y, 1);
if (nrowy ~= nobs)
    error('number of observations in y (%d) not equal to the number of rows of x (%d)',nrowy,nobs);
end

ne = options.dfmax;
if isempty(ne)
    ne = nvars + 1;
end

nx = options.pmax;
if isempty(nx)
    nx = min(ne * 2 + 20, nvars);
end

exclude = options.exclude;
if ~isempty(exclude)
    exclude = unique(exclude);
    if ~all(exclude > 0 & exclude <= nvars)
        error('Some excluded variables out of range');
    end
    jd = [length(exclude); exclude];
else
    jd = 0;
end
vp = options.penalty_factor;
if isempty(vp)
    vp = ones(1,nvars);
end

inparms = glmnetControl();

cl = options.cl;
if any(cl(1,:) > 0)
    error ('The lower bound should be non-positive');
end
if any(cl(2,:) < 0)
    error ('The lower bound should be non-negative');
end
cl(1,cl(1,:)==-Inf) = -inparms.big;
cl(2,cl(2,:)==Inf) = inparms.big;
if (size(cl,2) < nvars)
    if (size(cl,2) == 1)
        cl = cl * ones(1,nvars);
    else
        error('Require length 1 or nvars lower and upper limits');
    end
else
    cl = cl(:,1:nvars);
end

exit_rec = 0;
if (any(cl(:)==0))
    fdev = inparms.fdev;
    if (fdev ~= 0)
        optset.fdev = 0;
        glmnetControl(optset);
        exit_rec = 1;
    end
end

isd = double(options.standardize);
intr = double(options.intr);
if (intr == true) && (strcmp(family, 'cox'))
    warning('Cox model has no intercept');
end
jsd = options.standardize_resp;
thresh = options.thresh;
lambda = options.lambda;
lambda_min = options.lambda_min;

if isempty(lambda_min)
    if nobs < nvars
        lambda_min = 0.01;   
    else
        lambda_min = 1e-4;
    end
end

lempty = isempty(lambda);
if lempty
    if (lambda_min >= 1)
        error('lambda_min should be less than 1');
    end
    flmin = lambda_min;
    ulam = 0.0;
else
    flmin = 1.0;
    if any(lambda < 0)
        error ('lambdas should be non-negative');
    end
    ulam = sort(lambda,'descend');
    nlam = length(lambda);
end

maxit = options.maxit;

gtype = options.gtype; 
if isempty(gtype)
    if (nvars < 500)
        gtype = 'covariance';
    else
        gtype = 'naive';
    end
end
ltype = options.ltype;

indl = find(strncmp(ltype,{'Newton','modified.Newton'},length(ltype)),1);
if (isempty(indl))
    error('ltype should be one of ''Newton'', ''modified.Newton''');
else
    kopt = indl - 1;
end

if strcmp(family,'multinomial')
    mtype = options.mtype;
    indm = find(strncmp(mtype,{'ungrouped','grouped'},length(mtype)),1);
    if (isempty(indm))
        error('mtype should be one of ''ungrouped'', ''grouped''');
    else
        if (indm == 2)
            kopt = 2;
        end
    end
end

offset = options.offset;

is_sparse = false;
if issparse(x)
    is_sparse = true;
    [irs, jcs, x] = find(x);
    pcs = [0;cumsum(histc(jcs, 1:nvars))] + 1;
else
    irs = []; pcs = [];
end

if issparse(y)
    y = full(y);
end
  
switch family

    case 'gaussian'
        fit = elnet(x,is_sparse,irs,pcs,y,weights,offset,gtype,parm,lempty,...
            nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,family);
    case {'binomial', 'multinomial'}
        fit = lognet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
            jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,kopt,family);
    case 'cox'
        fit = coxnet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
            jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,maxit,family);
    case 'mgaussian'
        fit = mrelnet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
            jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,jsd,intr,maxit,family);
    case 'poisson'
        fit = fishnet(x,is_sparse,irs,pcs,y,weights,offset,parm,nobs,nvars,...
            jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,maxit,family);     
end

fit.call = {out_x, out_y, out_family, out_options};

if (exit_rec == 1)
    optset.fdev = fdev;
    glmnetControl(optset);
end

%------------------------------------------------------------------
% End function glmnet
%------------------------------------------------------------------