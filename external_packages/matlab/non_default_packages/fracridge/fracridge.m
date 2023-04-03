function [coef,alphas,offset] = fracridge(X,fracs,y,tol,mode,standardizemode)

% function [coef,alphas,offset] = fracridge(X,fracs,y,tol,mode,standardizemode)
%
% <X> is the design matrix (d x p) with different data points
%   in the rows and different regressors in the columns.
% <fracs> is a vector (1 x f) of one or more fractions between 0 and 1.
%   Fractions can be exactly 0 or exactly 1. However, values in between
%   0 and 1 should be no less than 0.001 and no greater than 0.999.
%   For example, <fracs> could be 0:.05:1 or 0:.1:1.
% <y> is the data (d x t) with one or more target variables in the columns.
% <tol> (optional) is a tolerance value such that singular values
%   below the tolerance are treated as 0. Default: 1e-6.
% <mode> (optional) can be:
%   0 means the default behavior
%   1 means to interpret <fracs> as exact alpha values that are desired.
%     This case involves using the rotation trick to speed up the ridge
%     regression but does not involve the fraction-based tailoring method.
%     In the case that <mode> is 1, the output <alphas> is returned as [].
%   Default: 0.
% <standardizemode> (optional) can be:
%   0 means the default behavior (do not modify any of the regressors). In this case, 
%     we do not add an offset term to the model. Note that the user may choose to include 
%     an constant regressor in <X> if so desired.
%   1 means to add an offset term to the model. In this case, the offset term is fully
%     estimated using ordinary least squares, and ridge regression is applied to 
%     de-meaned data and de-meaned regressors. (The user should not include a constant 
%     regressor in <X> if using this mode.)
%   2 means to standardize the regressors before performing ridge regression. In this case,
%     an offset term is added to the model, and is fully estimated using ordinary least
%     squares. Ridge regression is applied to de-meaned data and standardized
%     regressors. The returned regression weights will refer to the original regressors
%     (i.e., they will be adjusted for the effects of standardization). This mode may be
%     preferred for most applications given that the effects of regularization are 
%     influenced by the scale of the regressors. (The user should not include a constant
%     regressor in <X> if using this mode.)
%   Default: 0.
%
% return:
%  <coef> as the estimated regression weights (p x f x t)
%    for all fractional levels for all target variables.
%  <alphas> as the alpha values (f x t) that correspond to the
%    requested fractions. Note that alpha values can be Inf
%    in the case that the requested fraction is 0, and can
%    be 0 in the case that the requested fraction is 1.
%  <offset> as the offset term for each target variable (f x t).
%    Note that when <standardizemode> is 0, no offset term is added,
%    and <offset> is returned as all zeros.
%
% The basic idea is that we want ridge-regression solutions
% whose vector lengths are controlled by the user. The vector
% lengths are specified in terms of fractions of the length of
% the full-length solution. The full-length solution is either
% the ordinary least squares (OLS) solution in the case of
% full-rank design matrices, or the pseudoinverse solution in
% the case of rank-deficient design matrices.
%
% Framing the problem in this way provides several benefits:
% (1) we don't have to spend time figuring out appropriate
% alphas for each regression problem, (2) when implemented
% appropriately, we can compute the full set of solutions
% very efficiently and in a way that avoids the need
% to compute X'*X (which might be large), and (3) parameterizing
% the ridge-regression problem in terms of fractional lengths
% provides nicely interpretable results.
%
% The fraction-based method is tailored to each target variable
% in the sense that different target variables may need different
% alpha values in order to achieve the same fractional length.
% The output <alpha> is provided in case the user wants to know
% exactly what alpha values were used for each target variable.
%
% Notes:
% - We silently ignore regressors that are all zeros.
% - The class of <coef> and <alphas> is matched to the class of <X>.
% - All values in <X> and <y> should be finite. (A check for this is performed.)
% - If a given target variable consists of all zeros, this is
%   a degenerate case, and we will return regression weights that
%   are all zeros and alpha values that are all zeros.
%
% History:
% - 2020/07/26 - fix in interpolation to ensure monotonically increasing x-coordinates.
%                (code would have failed assert, prior to this.)
%
% % Example 1 (Demonstrate that fracridge achieves the correct fractional length)
%
% y = randn(1000,1);
% X = randn(1000,10);
% coef = inv(X'*X)*X'*y;
% [coef2,alpha] = fracridge(X,0.3,y);
% coef3 = inv(X'*X + alpha*eye(size(X,2)))*X'*y;
% coef4 = fracridge(X,alpha,y,[],1);
% norm(coef)
% norm(coef2)
% norm(coef2) ./ norm(coef)
% norm(coef2-coef3)
% norm(coef4-coef3)
%
% % Example 2 (Compare execution time between naive ridge regression and fracridge)
%
% y = randn(1000,300);        % 1000 data points x 300 target variables
% X = randn(1000,3000);       % 1000 data points x 3000 predictors
%   % naive approach
% tic;
% alphas = 10.^(-4:.5:5.5);  % guess 20 alphas
% cache1 = X'*X;
% cache2 = X'*y;
% coef = zeros(3000,length(alphas),300);
% for j=1:length(alphas)
%   coef(:,j,:) = permute(inv(cache1 + alphas(j)*eye(size(X,2)))*cache2,[1 3 2]);
% end
% toc;
%   % fracridge approach
% tic;
% fracs = .05:.05:1;         % get 20 equally-spaced lengths
% coef2 = fracridge(X,fracs,y);
% toc;
%   % fracridge implementation of simple rotation
% tic;
% coef3 = fracridge(X,alphas,y,[],1);
% toc;
% assert(all(abs(coef(:)-coef3(:))<1e-4));
%
% % Example 3 (Plot coefficient paths and vector length for a simple example)
%
% y = randn(100,1);
% X = randn(100,6)*(1+rand(6,6));
% fracs = .05:.05:1;
% [coef,alphas] = fracridge(X,fracs,y);
% figure;
% subplot(1,2,1); hold on;
% plot(fracs,coef');
% xlabel('Fraction');
% title('Trace plot of coefficients');
% subplot(1,2,2); hold on;
% plot(fracs,sqrt(sum(coef.^2,1)),'ro-');
% xlabel('Fraction');
% ylabel('Vector length');
%
% % Example 4 (Demonstrate how fracridge handles standardization of regressors)
%
% X = 20 + randn(50,2);
% y = X*rand(2,1) + randn(50,1);
% fracs = 0.1:0.1:1;
% [coef,alphas,offset] = fracridge(X,fracs,y,[],[],2);
% modelfit = X*coef + repmat(offset',[50 1]);
% figure; hold on;
% cmap0 = copper(length(fracs));
% h = [];
% legendstr = {};
% for p=1:length(fracs)
%   h(p) = plot(modelfit(:,p),'-','Color',cmap0(p,:));
%   legendstr{p} = sprintf('Frac %.1f',fracs(p));
% end
% hdata = plot(y,'k-','LineWidth',2);
% legend([hdata h],['Data' legendstr],'Location','EastOutside');

%% %%%%%% SETUP

% deal with inputs
if ~exist('tol','var') || isempty(tol)
  tol = 1e-6;
end
if ~exist('mode','var') || isempty(mode)
  mode = 0;
end
if ~exist('standardizemode','var') || isempty(standardizemode)
  standardizemode = 0;
end

% internal constants
bbig   = 10^3;    % sets huge-bias alpha. this will reduce length to less than 1/(1+10^3) = 0.000999 of the full length.
bsmall = 10^-3;   % sets small-bias alpha. this will preserve length to at least 1/(1+10^-3) = 0.999 of the full length.
bstep  = 0.2;     % step size for alpha in log10 space
debugmode = 0;    % if set to 1, we do some optional sanity checks

% sanity checks
assert(size(X,1)==size(y,1),'number of rows in X and y should match');
assert(ismember(class(X),{'single' 'double'}),'X should be in single or double format');
assert(ismember(class(y),{'single' 'double'}),'y should be in single or double format');
assert(all(isfinite(X(:))),'X must have only finite values');
assert(all(isfinite(y(:))),'y must have only finite values');

% ignore bad regressors (those that are all zeros)
bad = all(X==0,1);
if any(bad)
  fprintf('WARNING: bad regressors detected (all zero). ignoring these regressors.\n');
  X = X(:,~bad);
end

% deal with standardization
switch standardizemode

case 0

  % do nothing

case 1

  % de-mean the data and the regressors
  mnX = mean(X,1);
  mny = mean(y,1);
  X = X - repmat(mnX,[size(X,1) 1]);
  y = y - repmat(mny,[size(y,1) 1]);
  assert(~any(all(X==0,1)),'in <standardizemode>==1, there should not be a constant term in <X>');

case 2

  % standardize predictors and de-mean the data
  mnX = mean(X,1);
  sdX = std(X,[],1);
  mny = mean(y,1);
  X = (X - repmat(mnX,[size(X,1) 1])) ./ repmat(sdX,[size(X,1) 1]);
  y = y - repmat(mny,[size(y,1) 1]);
  assert(~any(sdX==0),'in <standardizemode>==2, there should not be a constant term in <X>');

end

% calc
d = size(X,1);      % number of data points
p = size(X,2);      % number of parameters (regressors)
t = size(y,2);      % number of distinct outputs (targets) being modeled
f = length(fracs);  % number of fractions (or alphas) being requested

%% %%%%%% PERFORM SVD AND ROTATE DATA

% decompose X [this is a costly step]
  % [u,s,v] = svd(X,'econ');  % when d>=p, u is d x p, s is p x p, v is p x p
  %                           % when d< p, u is d x d, s is d x d, v is p x d
if d >= p

  % avoid making a large u
  [~,s,v] = svd(X'*X,'econ');

  % extract the singular values
  selt = sqrt(diag(s));  % p x 1
  clear s;               % clean up to save memory

  % rotate the data, i.e. u'*y (X=u*s*v', u=X*v*inv(s), u'=inv(s)*v'*X')
  % in the below, the operation is (pxp) x (pxp) x (pxd) x (dxt).
  % we group appropriately to speed things up.
  if t >= d
    ynew = (diag(1./selt) * v' * X') * y;  % p x t
  else
    ynew = diag(1./selt) * v' * (X' * y);  % p x t
  end

else

  % do it
  [u,s,v] = svd(X,'econ');

  % extract the singular values
  selt = diag(s);  % d x 1
  clear s;         % clean up to save memory

  % rotate the data
  ynew = u'*y;     % d x t
  clear u;         % clean up to save memory

end

% calc
sz = length(selt);  % the size (rank) of the problem (either p or d)

% mark singular values that are essentially zero
isbad = selt < tol;      % p x 1 (OR d x 1)
anyisbad = any(isbad);
if anyisbad
  fprintf('WARNING: some singular values are being treated as 0.\n');
end

%% %%%%%% COMPUTE OLS SOLUTION IN ROTATED SPACE

% compute the OLS (or pseudoinverse) solution in the rotated space
ynew = ynew ./ repmat(selt,[1 t]);  % p x t (OR d x t)
if anyisbad
  ynew(isbad,:) = 0;  % the solution will be 0 along directions associated with singular values that are essentially zero
end

%% %%%%%% DO THE MAIN STUFF

% initialize
if d >= p
  coef = zeros(p,f*t,class(X));  % this is the easy case. the final rotation doesn't change the dimensionality.
else
  coef = zeros(d,f*t,class(X));  % in this case, the final rotation will change from d dimensions to p dimensions.
end
offset = zeros(f,t,class(X));

% we have two modes of operation...
switch mode

% this is the case of fractions being requested
case 0

  %% %%%%% DO SOME SETUP

  % figure out a reasonable grid for alpha at reasonable level of granularity
  val1 = bbig*selt(1)^2;              % huge bias (take the biggest singular value down massively)
  val2 = bsmall*min(selt(~isbad))^2;  % tiny bias (just add a small amount to the smallest singular value)
  alphagrid = fliplr([0 10.^(floor(log10(val2)):bstep:ceil(log10(val1)))]);  % huge bias to tiny bias to no bias (1 x g)
  g = length(alphagrid);

  % note that alphagrid is like [10e6 ... 0].
  %
  % also, note that the <bstep> is a heuristic. we just need to sample alphas at sufficient
  % granularity such that linear interpolation will be able to find the requested
  % fractions with sufficient accuracy. the idea here is that the computations for
  % constructing solutions with different regularization amounts (see below) are cheap,
  % and so we are happy to "oversample" to some extent.

  % construct scaling factor
  seltSQ = selt.^2;                                                     % p x 1 (OR d x 1)
  scLG = repmat(seltSQ,[1 g]);                                          % p x g (OR d x g)
  scLG = scLG ./ (scLG + repmat(alphagrid,[size(scLG,1) 1]));
  if anyisbad
    scLG(isbad,:) = 0;                                                  % for safety, ensure bad singular values get scalings of 0
  end

  % pre-compute for speed
  scLG = scLG.^2';                                                      % g x p (OR g x d)
  seltSQ2 = repmat(seltSQ,[1 f]);                                       % p x f (OR d x f)
  fracisz = find(fracs==0);                                             % indices into <fracs>
  logalpha = log(1+alphagrid)';                                         % transform alphas to log(1+x) scale (g x 1)

  % init
  alphas = zeros(f,t,class(X));

  %% %%%%% PROCEED TO COSTLY FOR-LOOP

  % compute ridge regression solutions and corresponding alphas
  for ii=1:t

    % compute vector length for each alpha in the grid
    len = sqrt(scLG*ynew(:,ii).^2);  % g x 1

    % when ynew is all 0, this is a degenerate case
    % and will cause len to be a bunch of zeros.
    % let's check the last element (no regularization),
    % and if it's 0, we know this is the degenerate case.
    % if so, we just continue, and this will result in alphas
    % and coeff to be left at 0.
    if len(end)==0
      continue;
    end

    % express lengths as fractions relative to the full length (i.e. no regularization)
    len = len / len(end);

    % inspection (is the proposed interpolation scheme reasonable?)
    if debugmode
      figure;
      scatter(len,logalpha,'ro');
      xlabel('length');
      ylabel('log(1+alphagrid)');
    end

    % sanity check that the gridding is not too coarse
    mxgap = max(abs(diff(len)));                    % maximum gap
    assert(mxgap < 0.2,'need to decrease bstep!');  % if this fails, bstep should be smaller

    % use linear interpolation to determine alphas that achieve the desired fractional levels.
    % we interpolate in log(1+x) space in order to help achieve good quality interpolation.
    iix = flipud([true; diff(flipud(len)) < -1e-5]);  % we need to ensure monotonically increasing x-coordinates!!
                                                      % note that we flipud to ensure that 1 is an "anchor" for this culling
    temp = interp1qr(len(iix),logalpha(iix),fracs');  % f x 1 (if out of range, will be NaN; we check this later)
    temp = exp(temp)-1;                             % undo the log transform
    temp(fracisz) = Inf;                            % when frac is exactly 0, we are out of range, so handle explicitly
    alphas(:,ii) = temp;

    % apply scaling to the OLS solution
    coef(:,(ii-1)*f+(1:f)) = repmat(seltSQ .* ynew(:,ii),[1 f]) ./ (seltSQ2 + repmat(temp',[sz 1]));

  end

  % accuracy check to see if the achieved vector lengths are close to what was requested
  if debugmode
    temp = sqrt(sum(coef.^2,1));    % 1 x f*t
    temp2 = sqrt(sum(ynew.^2,1));   % 1 x t
    temp3 = reshape(temp,[f t]) ./ repmat(temp2,[f 1]);  % f x t
    temp4 = repmat(fracs',[1 t]);
    figure; hold on;
    scatter(temp3(:),temp4(:),'r.');
    xlabel('empirical fraction');
    ylabel('requested fraction');
    plot([0 1],[0 1],'g-');
  end

  % for safety, ensure bad singular values get scalings of 0
  if anyisbad
    coef(isbad,:) = 0;
  end

  % if all went well, no value should be NaN (but can be Inf)
  assert(~any(isnan(alphas(:))),'NaN encountered in alphas. Is an element in <fracs> too close to 0?');


% this is the case of conventional alphas being requested
case 1

  % construct scaling factor
  sc = repmat(selt.^2,[1 f]);                                 % p x f (OR d x f)
  sc = sc ./ (sc + repmat(fracs,[size(sc,1) 1]));             % p x f (OR d x f)
  if anyisbad
    sc(isbad,:) = 0;                                          % for safety, ensure bad singular values get scalings of 0
  end

  % apply scaling to the OLS solutions.
  % do it in a for-loop to save memory usage.
  for ii=1:t
    coef(:,(ii-1)*f+(1:f)) = sc .* repmat(ynew(:,ii),[1 f]);
  end

  % deal with output (alphas is irrelevant, so set to [])
  alphas = cast([],class(X));

end

% rotate solution to the original space
if standardizemode==2
  coef = reshape((v./repmat(sdX',[1 size(v,2)]))*coef,[p f t]);  % adjust the regression weights for standardization effects on-the-fly
else
  coef = reshape(v*coef,[p f t]);
end

% calculate offset term
if ismember(standardizemode,[1 2])
  for ff=1:f
    offset(ff,:) = mny - mnX*permute(coef(:,ff,:),[1 3 2]);  % 1 x t - (1 x p)*(p x t)
  end
end

% deal with bad regressors
if any(bad)
  coef(~bad,:,:) = coef;
  coef(bad,:,:) = 0;
end
