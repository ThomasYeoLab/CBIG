function glmnetPlot( x, xvar, label, type, varargin )

%--------------------------------------------------------------------------
% glmnetPlot.m: plot coefficients from a "glmnet" object
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Produces a coefficient profile plot fo the coefficient paths for a
%    fitted "glmnet" object.
%
% USAGE:
%    glmnetPlot(fit);
%    glmnetPlot(fit, xvar);
%    glmnetPlot(fit, xvar, label);
%    glmnetPlot(fit, xvar, label, type);
%    glmnetPlot(fit, xvar, label, type, ...);
%    (Use empty matrix [] to apply the default value, eg. glmnetPlot(fit,
%    [], [], type).)
%
% INPUT ARGUMENTS:
% x           fitted "glmnet" model.
% xvar        What is on the X-axis. 'norm' plots against the L1-norm of
%             the coefficients, 'lambda' against the log-lambda sequence,
%             and 'dev' against the percent deviance explained.
% label       If true, label the curves with variable sequence numbers.
% type        If type='2norm' then a single curve per variable, else
%             if type='coef', a coefficient plot per response.
% varargin    Other graphical parameters to plot.
%
% DETAILS:
%    A coefficient profile plot is produced. If x is a multinomial model, a
%    coefficient plot is produced for each class.
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
%    glmnet, glmnetSet, glmnetPrint, glmnetPredict and glmnetCoef.
%
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    g2=randsample(2,100,true);
%    g4=randsample(4,100,true);
%    fit1=glmnet(x,y);
%    glmnetPlot(fit1);
%    glmnetPlot(fit1, 'lambda', true);
%    fit3=glmnet(x,g4,'multinomial');
%    glmnetPlot(fit3);
%
% DEVELOPMENT: 
%    14 Jul 2009: Original version of glmnet.m written.
%    30 Aug 2013: Updated glmnet.m with more options and more models
%                 (multi-response Gaussian, cox and Poisson models) supported.

if nargin < 2 || isempty(xvar)
    xvar = 'norm';
end

if nargin < 3 || isempty(label)
    label = false;
end

if nargin < 4 || isempty(type)
    type = 'coef';
end

xvarbase = {'norm','lambda','dev'};
xvarind = find(strncmp(xvar,xvarbase,length(xvar)),1);
if isempty(xvarind)
    error('xvar should be one of ''norm'', ''lambda'', ''dev''');
else
    xvar = xvarbase{xvarind};
end

typebase = {'coef','2norm'};
typeind = find(strncmp(type,typebase,length(type)),1);
if isempty(typeind)
    error('type should be one of ''coef'', ''2norm''');
else
    type = typebase{typeind};
end

if any(strcmp(x.class,{'elnet','lognet','coxnet','fishnet'}))
    plotCoef(x.beta,[],x.lambda,x.df,x.dev,label,xvar,'','Coefficients',varargin{:});
end

if strcmp(x.class,'multnet') || strcmp(x.class,'mrelnet')
    beta = x.beta;
    if strcmp(xvar,'norm')
        norm = 0;
        nzbeta = beta;
        for i=1:length(beta)
            which = nonzeroCoef(beta{i});
            nzbeta{i} = beta{i}(which,:);
            norm = norm + sum(abs(nzbeta{i}),1);
        end
    else
        norm = 0;
    end
    
    if strcmp(type,'coef')
        ncl = size(x.dfmat,1);
        if strcmp(x.class,'multnet')
            for i = 1:ncl
                plotCoef(beta{i},norm,x.lambda,x.dfmat(i,:),x.dev,label,xvar,'',sprintf('Coefficients: Class %d', i),varargin{:});
            end
        else
            for i = 1:ncl
                plotCoef(beta{i},norm,x.lambda,x.dfmat(i,:),x.dev,label,xvar,'',sprintf('Coefficients: Response %d', i),varargin{:});
            end
        end
    else
        dfseq = round(mean(x.dfmat,1));
        coefnorm = beta{1}*0;
        for i=1:length(beta)
            coefnorm = coefnorm + abs(beta{i}).^2;
        end
        coefnorm = sqrt(coefnorm);
        if strcmp(x.class,'multnet')
            plotCoef(coefnorm,norm,x.lambda,dfseq,x.dev,label,xvar,'',sprintf('Coefficient 2Norms'),varargin{:});
        else
            plotCoef(coefnorm,norm,x.lambda,x.dfmat(1,:),x.dev,label,xvar,'',sprintf('Coefficient 2Norms'),varargin{:});
        end
        
    end
end

%----------------------------------------------------------------
% End function glmnetPlot
%----------------------------------------------------------------

function plotCoef(beta,norm,lambda,df,dev,label,xvar,xlab,ylab,varargin)

which = nonzeroCoef(beta);
idwhich = find(which);  %row indices
nwhich = length(idwhich);
if nwhich == 0
    error('No plot produced since all coefficients zero')
end
if nwhich == 1
    warning('1 or less nonzero coefficients; glmnet plot is not meaningful');
end

beta = beta(which,:);
if strcmp(xvar, 'norm')
    if isempty(norm)
        index = sum(abs(beta),1);
    else
        index = norm;
    end
    iname = 'L1 Norm';
elseif strcmp(xvar, 'lambda')
    index=log(lambda);
    iname='Log Lambda';
elseif strcmp(xvar, 'dev')
    index=dev;
    iname='Fraction Deviance Explained';
end

if isempty(xlab)
    xlab = iname;
end

%Prepare for the figure (especially for the df labels)

clf;

beta = transpose(beta);
plot(index,beta,varargin{:});

axes1 = gca;
axes;
axes2 = gca;

xlim1 = get(axes1,'XLim');
ylim1 = get(axes1,'YLim');

%idxrange = range(index);
%atdf = linspace(min(index)+idxrange/12, max(index)-idxrange/12, 6);
atdf = get(axes1,'XTick');
indat = ones(size(atdf));
if (index(end) >= index(1))
    for j = length(index):-1:1
        indat(atdf <= index(j)) = j;
    end
else
    for j = 1:length(index)
        indat(atdf <= index(j)) = j;
    end
end
prettydf = df(indat);
prettydf(end) = df(end);

set(axes1,'box','off','XAxisLocation','bottom','YAxisLocation','left');
set(axes2,'XAxisLocation','top','YAxisLocation','right',...
    'XLim',[min(index),max(index)],'XTick',atdf,'XTickLabel',prettydf,...
    'YTick',[],'YTickLabel',[],'TickDir','out');
xlabel(axes2,'Degrees of Freedom')
axes(axes1);

line(xlim1,[ylim1(2),ylim1(2)],'Color','k');
line([xlim1(2),xlim1(2)],ylim1,'Color','k');

xlabel(xlab);
ylabel(ylab);


if (label)
    xpos = max(index);
    adjpos = 2;
    if strcmp(xvar,'lambda')
        xpos = min(index);
        adjpos = 1;
    end
    bsize = size(beta);
    for i = 1: bsize(2)
        text(1/2*xpos+1/2*xlim1(adjpos),beta(bsize(1),i),num2str(idwhich(i)));
    end
end

linkaxes([axes1 axes2],'x');

%----------------------------------------------------------------
% End private function plotCoef
%----------------------------------------------------------------
