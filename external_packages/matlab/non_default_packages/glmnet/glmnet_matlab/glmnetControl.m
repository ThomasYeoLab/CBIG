function ivals = glmnetControl(pars)

% DESCRIPTION:
%    View and/or change the factory default parameters in glmnet
%
% USAGE:
%   glmnetControl; (with no input or output arguments)
%   displays all inner parameters and their current values.
%   glmnetControl(pars);
%   sets the internal parameters that appear in the fields of pars to the
%   new values.
%
% ARGUMENTS:
% pars is a structure with the following fields.
% fdev        minimum fractional change in deviance for stopping path;
%             factory default = 1.0e-5.
% devmax      maximum fraction of explained deviance for stopping path;
%             factory default = 0.999.
% eps         minimum value of lambda.min.ratio (see glmnet); factory
%             default= 1.0e-6.
% big         large floating point number; factory default = 9.9e35. Inf in
%             definition of upper.limit is set to big.
% mnlam       minimum number of path points (lambda values) allowed;
%             factory default = 5.
% pmin        minimum null probability for any class. factory default =
%             1.0e-5.
% exmx        maximum allowed exponent. factory default = 250.0.
% prec        convergence threshold for multi response bounds adjustment
%             solution. factory default = 1.0e-10.
% mxit  	  maximum iterations for multiresponse bounds adjustment
%             solution. factory default = 100.
% factory     If true, reset all the parameters to the factory default;
%             default is false.
%
% DETAILS:
%    If called with no arguments, glmnetControl() returns a structure with 
%    the current settings of these parameters. Any arguments included in the
%    fields of the input structure sets those parameters to the new values, 
%    and then silently returns. The values set are persistent for the 
%    duration of the Matlab session.
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
%    glmnet.
%
% EXAMPLES:
%    pars = struct('fdev',0);
%    glmnetControl(pars);  %continue along path even though not much changes
%    glmnetControl();  %view current settings
%    pars = struct('factory',true);
%    glmnetControl(pars);  %reset all the parameters to their default

    if nargin == 0 || isempty(pars)
        [ivals.fdev, ivals.devmax, ivals.eps, ivals.big, ivals.mnlam, ...
            ivals.pmin, ivals.exmx, ivals.prec, ivals.mxit] = glmnetMex();
        if nargout == 0
            disp('internal paramters:');           
            disp( ivals );
        end
        return
    end
    
    if isfield(pars, 'factory') && (pars.factory == true)
        ivals.fdev = 1E-5;
        ivals.devmax = 0.999;
        ivals.eps = 1E-6;
        ivals.big = 9.9E+35;
        ivals.mnlam = 5;
        ivals.pmin = 1E-5;
        ivals.exmx = 250;
        ivals.prec = 1E-10;
        ivals.mxit = 100;
        task = 0;
        glmnetMex(task, ivals.fdev, ivals.devmax, ivals.eps, ivals.big, ivals.mnlam, ...
        ivals.pmin, ivals.exmx, ivals.prec, ivals.mxit);
        
    else
        [ivals.fdev, ivals.devmax, ivals.eps, ivals.big, ivals.mnlam, ...
            ivals.pmin, ivals.exmx, ivals.prec, ivals.mxit] = glmnetMex();
        vfields = fieldnames(ivals);
        for i = 1:(length(vfields))
            field = vfields{i};
            if isfield(pars, field)
                ivals.(field) = pars.(field);
            end
        end
        disp( ivals );
        task = 0;
        glmnetMex(task, ivals.fdev, ivals.devmax, ivals.eps, ivals.big, ivals.mnlam, ...
        ivals.pmin, ivals.exmx, ivals.prec, ivals.mxit);
    end
end