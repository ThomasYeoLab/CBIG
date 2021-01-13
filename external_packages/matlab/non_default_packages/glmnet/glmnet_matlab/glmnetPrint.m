function glmnetPrint( x )

%--------------------------------------------------------------------------
% glmnetPrint.m: print a glmnet object
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Print a summary of the glmnet path at each step along the path.
%
% USAGE: 
%    glmnetPrint(fit)
%
% INPUT ARGUMENTS:
% x         fitted glmnet object
%
% DETAILS:
%    Three-column matrix with columns Df, %Dev and Lambda is printed. The Df
%    column is the number of nonzero coefficients (Df is a reasonable name
%    only for lasso fits). %Dev is the percent deviance explained (relative
%    to the null deviance).
%
% LICENSE: GPL-2
%
% DATE: 14 Jul 2009
%
% AUTHORS:
%    Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani
%    Fortran code was written by Jerome Friedman
%    R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
%    The original MATLAB wrapper was written by Hui Jiang (14 Jul 2009),
%    and was updated and maintained by Junyang Qian (30 Aug 2013) junyangq@stanford.edu,
%    Department of Statistics, Stanford University, Stanford, California, USA.
%
% REFERENCES:
%    Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, 
%    http://www.jstatsoft.org/v33/i01/
%    Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
%
% SEE ALSO:
%    glmnet, glmnetSet, glmnetPredict and glmnetCoef methods.
% 
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    fit1=glmnet(x,y);
%    glmnetPrint(fit1);
%
% DEVELOPMENT: 
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported.

disp(sprintf('\tDf\t%%Dev\t\tLambda'));
% disp([x.df, x.dev, x.lambda]);
for i=1:length(x.lambda)
    disp(sprintf('%d\t%d\t%f\t%f', i, x.df(i), x.dev(i), x.lambda(i)));
end
