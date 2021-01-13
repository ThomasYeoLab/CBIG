function CVpred = cvglmnetPredict(object, newx, s, varargin)

%--------------------------------------------------------------------------
% cvglmnetPredict makes predictions from a "cv.glmnet" object.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    This function makes predictions from a cross-validated glmnet model,
%    using the stored "glmnet_fit" object, and the optimal value chosen for
%    lambda.
%
% USAGE:
%    pred = cvglmnetPredict(cvfit)
%    pred = cvglmnetPredict(cvfit, newx)
%    pred = cvglmnetPredict(cvfit, newx, s)
%    pred = cvglmnetPredict(cvfit, newx, s, ...)
%
% INPUT ARGUMENTS:
% object      Fitted "glmnet" model object.
% newx        Matrix of new values for x at which predictions are to be
%             made. Must be a matrix; can be sparse. See documentation for 
%             glmnetPredict.
% s           Value(s) of the penalty parameter lambda at which predictions
%             are required. Default is the value s='lambda_1se' stored on
%             the CV object. Alternatively s='lambda_min' can be used. If s
%             is numeric, it is taken as the value(s) of lambda to be used.
% varargin    Other arguments to predict.
%
% OUTPUT ARGUMENTS:
%    If only the cv.glmnet is provided, the function returns the 
%    coefficients at the default s = 'lambda_1se'. Otherwise, the object 
%    returned depends the ... argument which is passed on to the 
%    glmnetPredict for glmnet objects.
%             
%
% DETAILS:
%    This function makes it easier to use the results of cross-validation
%    to make a prediction. 
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
%    cvglmnet and glmnetPredict.
%
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    cvfit=cvglmnet(x,y);
%    pred1 = cvglmnetPredict(cvfit,x(1:5,:));
%    pred2 = cvglmnetPredict(cvfit,x(1:5,:), [0.001;0.002]);
%
% DEVELOPMENT:
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported. 

% s is a numeric value or either 'lambda.1se' or 'lambda.min'

if nargin < 2
    CVpred = cvglmnetCoef(object);
    return;
end

if nargin < 3 || isempty(s)
    s = 'lambda_1se';
end

if isnumeric(s)
    lambda = s;
else
    if any(strcmp(s, {'lambda_1se','lambda_min'}))
        lambda = object.(s);
    else
        error('Invalid form for s');
    end
end

CVpred = glmnetPredict(object.glmnet_fit,newx,lambda,varargin{:});

end
        
    