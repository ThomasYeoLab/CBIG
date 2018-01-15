function [h,pValue,stat,cValue] = lratiotest(uLL,rLL,dof,alpha)

% [h,pValue,stat,cValue] = lratiotest(uLL,rLL,dof,alpha)
%LRATIOTEST Likelihood ratio test of model specification
%
% Syntax:
%
%   [h,pValue,stat,cValue] = lratiotest(uLL,rLL,dof)
%   [h,pValue,stat,cValue] = lratiotest(uLL,rLL,dof,alpha)
%
% Description:
%
%   The likelihood ratio test compares specifications of nested models by
%   assessing the significance of restrictions to an extended model with
%   unrestricted parameters. Loglikelihoods maximized with respect to
%   restricted and unrestricted model parameters (rLL and uLL,
%   respectively) are used to compute the test statistic
%
%     stat = 2*(uLL - rLL)
%
%   When the test statistic exceeds a critical value in its asymptotic
%   distribution, the test rejects the null, restricted model in favor of
%   the alternative, unrestricted model. The asymptotic distribution is
%   chi-square, with degree-of-freedom parameter (dof) equal to the number
%   of restrictions. The nominal significance level of the test (alpha)
%   determines the critical value.
%
% Input Arguments:
%
%   uLL - Loglikelihoods optimized with respect to parameters for the
%       unrestricted models to be tested. If uLL is a scalar, it is
%       expanded to the same length as rLL. If uLL and rLL are both
%       vectors, they must be the same length. If uLL is a row vector,
%       output arguments are also row vectors.
%
%   rLL - Loglikelihoods optimized with respect to parameters for the
%       restricted models to be tested. If rLL is a scalar, it is expanded
%       to the same length as uLL. If rLL and uLL are both vectors, they
%       must be the same length. Elements of rLL should not be greater than
%       corresponding elements of uLL. If rLL is a row vector, output
%       arguments are also row vectors. 
%
%   dof - Degree-of-freedom parameters for the asymptotic chi-square
%       distributions of the test statistics. Elements of dof are positive
%       integers equal to the number of restrictions in the corresponding
%       model comparison, and should be less than the number of parameters
%       in the unrestricted model. If dof is a scalar, it is expanded to a
%       vector with length equal to the number of tests. If dof is a
%       vector, it must have length equal to the number of tests.
%
% Optional Input Arguments:
%
%   alpha - Nominal significance levels for the tests. Elements of alpha
%       must be greater than zero and less than one. If alpha is a scalar,
%       it is expanded to a vector with length equal to the number of
%       tests. If alpha is a vector, it must have length equal to the
%       number of tests. The default value of alpha is 0.05.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to the
%       number of tests. Values of h equal to 1 indicate rejection of the
%       null, restricted model in favor of the alternative, unrestricted
%       model. Values of h equal to 0 indicate a failure to reject the
%       restricted model.
%
%   pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests.
%
%   stat - Vector of test statistics, with length equal to the number of
%       tests.
%
%   cValue - Vector of critical values for the tests, determined by alpha,
%       with length equal to the number of tests.
%
% Notes:
%
%   o LRATIOTEST performs multiple, independent tests when either uLL or
%     rLL is a vector. If uLL is a scalar and rLL is a vector, LRATIOTEST
%     "tests down" against multiple restricted models. If rLL is a scalar
%     and uLL is a vector, LRATIOTEST "tests up" against multiple
%     unrestricted models. If both uLL and rLL are vectors, LRATIOTEST
%     compares model specifications pairwise.
%   
%   o The functions ARIMA/ESTIMATE and ARIMA/INFER return loglikelihoods of
%     time series data using GARCH models. Use these outputs as inputs to
%     LRATIOTEST to compare specifications of GARCH(p,q) models.
%
%   o The significance level alpha of LRATIOTEST is nominal, in that it
%     specifies a rejection probability in the asymptotic distribution. The
%     actual rejection probability will generally be greater than the
%     nominal significance. 
%
%   o Likelihood ratio tests are useful when both unrestricted and
%     restricted parameter estimates are easily computed. By comparison,
%     the Wald test (WALDTEST) requires only unrestricted parameter
%     estimates and the Lagrange multiplier test (LMTEST) requires only
%     restricted parameter estimates.
%
% Example:
%
%   % Load exchange rate data:
%   load Data_MarkPound
%   returns = price2ret(Data);
% 
%   % Fit GARCH(1,1) and GARCH(2,1) models to the data:
%   spec1 = garch(1,1);
%   [~,~,rLL] = estimate(spec1,returns);
%   spec2 = garch(2,1);
%   [~,~,uLL] = estimate(spec2,returns);
% 
%   % Compare model specifications:
%   [h,pValue,stat,cValue] = lratiotest(uLL,rLL,1)
% 
%   % The test fails to reject spec1 in favor of spec2.
%
% References:
%
%   [1] Davidson, R. and J. G. MacKinnon. Econometric Theory and Methods.
%       Oxford, UK: Oxford University Press, 2004.
%   
%   [2] Godfrey, L. G. Misspecification Tests in Econometrics. Cambridge,
%       UK: Cambridge University Press, 1997.
% 
%   [3] Greene, W. H. Econometric Analysis. Upper Saddle River, NJ: Pearson
%       Prentice Hall, 2008.
% 
%   [4] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also WALDTEST, LMTEST, ARIMA/ESTIMATE, ARIMA/INFER.


% Copyright 2009-2010 The MathWorks, Inc.

% Check for sufficient inputs:

if (nargin < 3) || isempty(uLL) || isempty(rLL) || isempty(dof)
    
    error(message('econ:lratiotest:UnspecifiedInput'))
      
end

% Check uLL and rLL inputs:

if ~isvector(uLL) || ~isvector(rLL)
    
    error(message('econ:lratiotest:LLNonVector'))
      
else
    scalarExpansion = true; % Flag to expand scalar inputs
    rowOutput = false;      % Flag to display row outputs
    
    if isscalar(uLL) && ~isscalar(rLL)
        rowOutput = (size(rLL,1) == 1);
        if rowOutput
            rLL = rLL(:); % Column for vector comparisons
        end
        uLL = uLL(ones(size(rLL))); % Scalar expansion
        
    elseif ~isscalar(uLL) && isscalar(rLL)
        rowOutput = (size(uLL,1) == 1);
        if rowOutput
            uLL = uLL(:); % Column for vector comparisons
        end
        rLL = rLL(ones(size(uLL))); % Scalar expansion
                
    elseif ~isscalar(uLL) && ~isscalar(rLL)
        if length(uLL) ~= length(rLL)
            
            error(message('econ:lratiotest:LLLengthMismatch'))
              
        end
        
        rowULL = (size(uLL,1) == 1);
        rowRLL = (size(rLL,1) == 1);
        rowOutput = rowULL || rowRLL;
        
        if rowULL
            uLL = uLL(:); % Column for vector comparisons
        end
        
        if rowRLL
            rLL = rLL(:); % Column for vector comparisons
        end
        
    else % Both LLs are scalars
        scalarExpansion = false;
        
    end
    
end

numTests = length(uLL); % Number of tests

if any(rLL > uLL)
    
    warning(message('econ:lratiotest:RLLExceedsULL'))
        
end

% Check dof input:

if ~isvector(dof)
    
    error(message('econ:lratiotest:DofNonVector'))
      
elseif any(mod(dof,1) ~= 0) || any(dof <= 0)
    
    error(message('econ:lratiotest:DofNonPositiveInteger'))
          
elseif ~isscalar(dof) && (length(dof) ~= numTests)
    
    error(message('econ:lratiotest:DofLengthMismatch'))
 
elseif isscalar(dof) && scalarExpansion
    
    dof = dof(ones(size(uLL))); % Scalar expansion
      
elseif ~isscalar(dof)
    
    dof = dof(:); % Column for vector comparisons
    
end

% Check alpha input or set default:

if (nargin >= 4) && ~isempty(alpha)
    
    if ~isvector(alpha)
        
        error(message('econ:lratiotest:AlphaNonVector'))
                    
    elseif any(alpha <= 0) || any(alpha >= 1)
        
        error(message('econ:lratiotest:AlphaOutOfRange'));
        
    elseif ~isscalar(alpha) && (length(alpha) ~= numTests)
        
        error(message('econ:lratiotest:AlphaLengthMismatch'));
          
    elseif isscalar(alpha) && scalarExpansion
        
        alpha = alpha(ones(size(uLL))); % Scalar expansion
          
    elseif ~isscalar(alpha)
        
        alpha = alpha(:); % Column for vector comparisons
        
    end
   
else % Set default
   alpha = 0.05;
   alpha = alpha(ones(size(uLL))); % Scalar expansion
   
end

% Perform the tests:

stat = 2*(uLL-rLL);
pValue = 1-chi2cdf(stat,dof);
h = (pValue <= alpha);

if nargout >= 4
    cValue = chi2inv(1-alpha,dof);
else
    cValue = [];
end

% Display outputs as row vectors if either uLL or rLL is a row vector:

if rowOutput
    h = h';
    pValue = pValue';
    stat = stat';
    cValue = cValue';
end