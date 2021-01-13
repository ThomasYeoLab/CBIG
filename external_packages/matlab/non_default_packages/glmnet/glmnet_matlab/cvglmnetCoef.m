function result = cvglmnetCoef(object, s)

%--------------------------------------------------------------------------
% cvglmnetCoef computes coefficients from a "cv.glmnet" object.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    This function extracts coefficients at certain lambdas if they are
%    in the lambda sequence of a "cv.glmnet" object or make predictions
%    if they are not.
%
% USAGE:
%    mcoef=cvglmnetCoef(object);
%    ncoef=cvglmnetCoef(object, s);
%
% INPUT ARGUMENTS:
% object      Fitted "glmnet" model object.
% s           Value(s) of the penalty parameter lambda at which computation
%             is required. Default is the value s='lambda_1se' stored on
%             the CV object. Alternatively s='lambda_min' can be used. If s
%             is numeric, it is taken as the value(s) of lambda to be used.
%
% OUTPUT ARGUMENTS:
% result      If s is 'lambda_1se' or 'lambda_min', the coefficients at 
%             that s is returned. If s is numeric, a (nvars+1) x length(s) 
%             matrix is returned with each column being the coefficients 
%             at an s. Note that the first row are the intercepts (0 if no 
%             intercept in the original model).
%
% DETAILS:
%    The function uses linear interpolation to make predictions for values 
%    of s that do not coincide with those used in the fitting algorithm. 
%    Exact prediction is not supported currently.
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
%    cvglmnet, cvglmnetPrint, and cvglmnetPredict.
%
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    cvfit=cvglmnet(x,y);
%    ncoef=cvglmnetCoef(cvfit,'lambda_min');
%
% DEVELOPMENT:
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported. 

if nargin < 2 || isempty(s)
    s = object.lambda_1se;
end

if isnumeric(s)
    lambda = s;
elseif ischar(s)
    sbase = {'lambda_1se','lambda_min'};
    sind = find(strncmp(s,sbase,length(s)),1);
    s = sbase{sind};
    lambda = object.(s);
else
    error('Invalid form of s');
end

result = glmnetCoef(object.glmnet_fit, lambda);

end