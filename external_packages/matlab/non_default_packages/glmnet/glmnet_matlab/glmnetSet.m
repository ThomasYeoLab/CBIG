function options = glmnetSet(opts)

%--------------------------------------------------------------------------
% glmnetSet creates or alters an options structure for glmnet.m.
%--------------------------------------------------------------------------
%   options = glmnetSet; (with no input arguments)
%   creates a structure with all fields set to their default values.
%   Each field is an option (also called a parameter).
%
%   glmnetSet; (with no input or output arguments)
%   displays all options and their default values.
%
%   options = glmnetSet(opts); 
%   creates a structure with all fields set to their default values,
%   except valid fields in the structure "opts" replace the defaults.
%
% options.alpha       The elasticnet mixing parameter, with 0 < alpha <= 1.
%                     The penalty is defined as
%                           (1-alpha)/2(||beta||_2)^2+alpha||beta||_1.
%                     Default is alpha = 1, which is the lasso penalty;
%                     Currently alpha = 0 the ridge penalty.
% options.nlambda     The number of lambda values - default is 100.
% options.lambda      A user supplied lambda sequence. Typical usage is to
%                     have the program compute its own lambda sequence
%                     based on nlambda and lambda_min. Supplying a value of
%                     lambda override this. WARNING: Use with care. Do not 
%                     supply a single value for lambda (for predictions 
%                     after CV use cvglmnetPredict() instead). Supply a 
%                     decreasing sequence of lambda values. glmnet relies
%                     on its warm starts for speed, and it's often faster
%                     to fit a whole path than compute a single fit.
% options.standardize Logical flag for x variable standardization, prior to
%                     fitting the model sequence. The coefficients are
%                     always returned on the original scale. Default is
%                     standardize = true. If variables are in the same
%                     units already, you might not wish to standardize. See
%                     details below for y standardization with
%                     family='gaussian'.
% options.weights     Observation weights. Can be total counts if responses
%                     are proportion matrices. Default is 1 for each
%                     observation.
% options.intr        Should intercept(s) be fitted (default=true) or set
%                     to zero (false).
% options.offset      A vector of length nobs that is included in the
%                     linear predictor (a nobs x nc matrix for the
%                     "multinomial" family). Useful for the "poisson"
%                     family (e.g. log of exposure time), or for refining a
%                     model by starting at a current fit. Default is []. If
%                     supplied, then values must also be supplied to the
%                     predict function.
% options.lambda_min  Smallest value for lambda, as a fraction of
%                     lambda_max, the (data derived) entry value (i.e., the
%                     smallest value for which all coefficients are zero).
%                     The default depends on the sample size nobs relative
%                     to the number of variables nvars. If nobs > nvars,
%                     the default is 0.0001, close to zero. If nobs <
%                     nvars, the defaults is 0.01. A very small value of
%                     lambda_min will lead to a saturated fit. This is
%                     undefined for "binomial" and "multinomial" models,
%                     and glmnet will exit gracefully when the percentage
%                     deviance explained is almost 1.
% options.thresh      Convergence threshold for coordinate descent. Each 
%                     inner coordinate-descent loop continues until the 
%                     maximum change in the objective after any coefficient 
%                     update is less than thresh times the null deviance. 
%                     Defaults value is 1E-4.
% options.dfmax       Limit the maximum number of variables in the model. 
%                     Useful for very large nvars, if a partial path is
%                     desired. Default is nvars + 1.
% options.pmax        Limit the maximum number of variables ever to be
%                     nonzero. Default is min(dfmax * 2 + 20, nvars).
% options.exclude     Indices of variables to be excluded from the model. 
%                     Default is none. Equivalent to an infinite penalty
%                     factor (next item).
% options.penalty_factor
%                     Separate penalty factors can be applied to each
%                     coefficient. This is a number that multiplies lambda
%                     to allow differential shrinkage. Can be 0 for some
%                     variables, which implies no shrinkage, and that
%                     variable is always included in the model. Default is
%                     1 for all variables (and implicitly infinity for
%                     variables listed in exclude). Note: the penalty
%                     factors are internally rescaled to sum to nvars, and
%                     the lambda sequence will reflect this change.
% options.maxit       Maximum number of passes over the data for all lambda
%                     values; default is 10^5.
% options.cl          Two-row matrix with the first row being the lower 
%                     limits for each coefficient and the second the upper
%                     limits. Can be presented as a single column (which
%                     will then be replicated), else a matrix of nvars
%                     columns. Default [-Inf;Inf].
% options.gtype       Two algorithm types are supported for (only)
%                     family = 'gaussian'. The default when nvar<500 is
%                     options.gtype = 'covariance', and saves all
%                     inner-products ever computed. This can be much faster
%                     than options.gtype='naive', which loops through nobs
%                     every time an inner-product is computed. The latter
%                     can be far more efficient for nvar >> nobs
%                     situations, or when nvar > 500.
% options.ltype       If 'Newton' then the exact hessian is used (default),
%                     while 'modified.Newton' uses an upper-bound on the
%                     hessian, and can be faster.
% options.standardize_resp
%                     This is for the family='mgaussian' family, and allows
%                     the user to standardize the response variables.
% options.mtype       If 'grouped' then a grouped lasso penalty is used on
%                     the multinomial coefficients for a variable. This
%                     ensures they are all in our out together. The default
%                     is 'ungrouped'.  
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
%    and was updated and maintained by Junyang Qian (30 Aug 2013) junyangq@stanford.edu,
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
%    glmnet, cvglmnet methods.
% 
% DEVELOPMENT: 
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported.

% Set default options.
  options.weights          =           [];
  options.offset           =           [];        
  options.alpha            =          1.0;
  options.nlambda          =          100;
  options.lambda_min       =           [];
  options.lambda           =           [];
  options.standardize      =         true;
  options.intr             =         true;
  options.thresh           =         1E-7;
  options.dfmax            =           [];
  options.pmax             =           [];
  options.exclude          =           [];
  options.penalty_factor   =           [];
  options.cl               =   [-Inf;Inf];
  options.maxit            =         1E+5;
  options.gtype            =           [];        
  options.ltype            =     'Newton';
  options.standardize_resp =        false;
  options.mtype =             'ungrouped'; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End default options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quick return if no user opts
  if nargin == 0 || isempty(opts)
    if nargout == 0    % Display options.
      disp('pdco default options:')
      disp( options )
    end
    return
  end

% List of valid field names
  vfields = fieldnames( options );

% Grab valid fields from user's opts
  for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opts, field );
      options.(field) = opts.(field);
    end
  end
