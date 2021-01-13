function cvglmnetPlot(cvobject,sign_lambda,varargin)
%--------------------------------------------------------------------------
% cvglmnetPlot.m: plot the cross-validation curve produced by cvglmnet
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Plots the cross-validation curve, and upper and lower standard
%    deviation curves, as a function of the lambda values used. 
%
% USAGE:
%    cvglmnetPlot(cvfit);
%    cvglmnetPlot(cvfit, sign_lambda);
%    cvglmnetPlot(cvfit, sign_lambda, varagin);
%    (Use empty matrix [] to apply the default value, eg.
%    cvglmnetPlot(cvfit, [], 'linewidth', 2)).
%
% INPUT ARGUMENTS:
% cvobject    fitted "cv.glmnet" object
% sign_lambda Either plot against log(lambda) (default) or its negative if
%             sign_lambda=-1. 
% varargin    Other errorbar parameters.
% 
% DETAILS:
%    A plot is produced, and nothing is returned.
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
% SEE ALSO:
%    cvglmnet and glmnet.
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
%
%    cvob2=cvglmnet(x,ly,'binomial');
%    cvglmnetPlot(cvob2);
%    figure;
%
%    cvob3=cvglmnet(x,ly,'binomial',[],'class');
%    cvglmnetPlot(cvob3);   
%
% DEVELOPMENT: 
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported. 

    if nargin < 2 || isempty(sign_lambda)
        sign_lambda = 1;
    end
    
    sloglam = sign_lambda*log(cvobject.lambda);
    
    figure;

    errorbar(sloglam,cvobject.cvm,cvobject.cvsd,'Color',[0.5 0.5 0.5],varargin{:});
    hold on
    plot(sloglam,cvobject.cvm,'LineStyle','-','Marker','o','Color','r');
    axes1 = gca;
    xlim1 = get(axes1,'XLim');
    ylim1 = get(axes1,'YLim');
    line(sign_lambda*log([cvobject.lambda_min cvobject.lambda_min]),ylim1,'Color','b','LineStyle','--',...
            'linewidth',1)

    if cvobject.lambda_min ~=cvobject.lambda_1se
        line(sign_lambda*log([cvobject.lambda_1se cvobject.lambda_1se]),ylim1,'Color','b','LineStyle','--',...
            'linewidth',1)
    end
    
    axes;
    axes2 = gca;
    
    atdf = linspace(min(sloglam), max(sloglam), 12);
    
    indat = ones(size(atdf));
    if (sloglam(end) >= sloglam(1))
        for j = length(sloglam):-1:1
            indat(atdf <= sloglam(j)) = j;
        end
    else
        for j = 1:length(sloglam)
            indat(atdf <= sloglam(j)) = j;
        end
    end
    
    prettydf = cvobject.nzero(indat);
    
    set(axes1,'box','off','XAxisLocation','bottom','YAxisLocation','left');
    set(axes2,'XAxisLocation','top','YAxisLocation','right',...
        'XLim',xlim1,'XTick',atdf,'XTickLabel',prettydf,...
        'YTick',[],'YTickLabel',[],'TickDir','out');
    xlabel(axes2,'Degrees of Freedom')
    axes(axes1);
    
    line(xlim1,[ylim1(2),ylim1(2)],'Color','k');
    line([xlim1(2),xlim1(2)],ylim1,'Color','k');
    
    if (sign_lambda < 0)
        xlabel('-log(Lambda)');
    else
        xlabel('log(Lambda)');
    end
    ylabel(cvobject.name);
    
    linkaxes([axes1 axes2],'x');

    hold off
end