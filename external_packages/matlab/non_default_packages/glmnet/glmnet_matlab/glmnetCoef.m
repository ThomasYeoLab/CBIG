function result = glmnetCoef(object, s, exact)

%--------------------------------------------------------------------------
% glmnetCoef computes coefficients from a "glmnet" object.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    This function extracts coefficients at certain lambdas if they are
%    in the lambda sequence of a "glmnet" object or make predictions
%    if they are not in that sequence.
%
% USAGE:
%    glmnetCoef(object, s, exact)
%
%    Fewer input arguments(more often) are allowed in the call, but must
%    come in the order listed above. To set default values on the way, use
%    empty matrix []. 
%    For example, ncoef = glmnetCoef(fit,[],false).
%
% INPUT ARGUMENTS:
% object      Fitted "glmnet" model object.
% s           Value(s) of the penalty parameter lambda at which computation
%             is required. Default is the entire sequence used to create
%             the model.
% exact       If exact=true, and computation is to be made at values of s 
%             not included in the original fit, these values of s are merged
%             with object.lambda, and the model is refit before predictions
%             are made. If exact=false (default), then the function uses
%             linear interpolation to make predictions for values of s
%             that do not coincide with those used in the fitting
%             algorithm. Note that exact=true is fragile when used inside a
%             nested sequence of function calls. glmnetCoef() needs to
%             update the model, and expects the data used to create it in
%             the parent environment.
%
% OUTPUT ARGUMENTS:
% result      A (nvars+1) x length(s) matrix with each column being the 
%             coefficients at an s. Note that the first row are the 
%             intercepts (0 if no intercept in the original model).
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
%    glmnet, glmnetPrint, glmnetPredict, and cvglmnet.
%
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    fit=glmnet(x,y);
%    ncoef=glmnetCoef(fit,0.01,true);
%
% DEVELOPMENT:
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported. 

if nargin < 2
    s = object.lambda;
end

if nargin < 3
    exact = false;
end

if (exact && ~isempty(s))
    which = ismember(s,object.lambda);
    if ~all(which)
        lambda = unique([object.lambda;reshape(s,length(s),1)]);
        %-----create a new variable in the parent environment
        vname = 'newlam';
        expr = sprintf('any(strcmp(''%s'', who))',vname);
        newname = vname;
        i = 0;
        while (evalin('caller',expr))
            i = i + 1;
            newname = [vname,num2str(i)];
            expr = sprintf('any(strcmp(who,''%s''))',newname);
        end
        parlam = newname;
        %-----
        assignin('caller', parlam, lambda);
        
        vname = 'temp_opt';
        expr = sprintf('any(strcmp(''%s'', who))',vname);
        newname = vname;
        i = 0;
        while (evalin('caller',expr))
            i = i + 1;
            newname = [vname,num2str(i)];
            expr = sprintf('any(strcmp(who,''%s''))',newname);
        end
        paropt = newname;
        
        if strcmp('[]',object.call{3})
            famcall = object.call{3};
        else
            famcall = sprintf('''%s''',object.call{3});
        end
        
        if ~strcmp('[]', object.call{4})
            evalin('caller', strcat(paropt,'=',object.call{4},';'));
            evalin('caller', strcat(paropt,'.lambda = ',parlam,';'));
            newcall = sprintf('glmnet(%s, %s, %s, %s)', ...
                object.call{1}, object.call{2}, famcall, paropt);
            object = evalin('caller', newcall);
        else
            evalin('caller', strcat(paropt,'.lambda = ',parlam,';'));
            newcall = sprintf('glmnet(%s, %s, %s, %s)', ...
                object.call{1}, object.call{2}, famcall, paropt);
            object = evalin('caller', newcall);
        end
        evalin('caller', sprintf('clearvars %s %s;',parlam,paropt));
    end
end

result = glmnetPredict(object,[],s,'coefficients');

end