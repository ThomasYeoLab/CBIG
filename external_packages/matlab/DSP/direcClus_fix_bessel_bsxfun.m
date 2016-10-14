
function clustered = direcClus_fix_bessel_bsxfun(x,k,d,noruns,lambda,p,mInit,epsilon,flagCorStop, max_iter, assignment_flag)

% clustered = DirClus(x,k,d,noruns,lambda,p,mInit,epsilon,flagCorStop)
% 
% NOTE: This is a modified version of Danial Lashkari's original code direcClus.m,
%       which can be found here:
%       http://people.csail.mit.edu/danial/Site/Code.html
%
%       Our modification includes:
%       - Fixing the problem for the inconsistency of Bessel function in old
%       and new MATLAB versions.
%       - Changing all places that using repmat() to bsxfun()
%       - Changing breaking conditions
%       - Slightly difference of variable names, initialization convention
%
% d:              the dimensionality. If enter zero, the actual dimensionality
%                 of the data.
% noruns:         number of repetitions when using random initializations.
% initial values: initial concentration parameter (lambda), initial weights (p),
%                 initial cluster centers (mInit). If enter zero,
%                 initialized by random values.
% epsilon:        value of the condition for breaking the the iterations
% flagCorStop:    type of breaking condition. if set 1, uses the correlation
%                 between the updated means and the old ones
%                 in each iteration as the condition for breaking the loop, and does so by
%                 having such correlations for all the k mean vectors to be at least 1-epsilon.
%
%

%--------------------------------------------------------------------------
%    Initialization
%--------------------------------------------------------------------------

clustered.listlik = [];
clustered.listavcor = [];

% Find the parameters

n = size(x,1);
if d==0
    d = size(x,2);
end
clustered.d = d;
alpha = d/2 - 1;

% Normalize the data to the sphere

x = bsxfun(@times, x, 1./sqrt(sum(x.^2,2)) );

% Check the input variables

if ~exist('flagCorStop')
    flagCorStop = 0;
end

if ~exist('noruns')
    noruns = 1;
end

mInitFlag = 1;
if ~exist('mInit')
    mInitFlag = 0;
else if mInit == 0;
        mInitFlag = 0;
    end
end

epsilonCheck = 0;
if exist('epsilon')
    epsilonCheck = 1;
end

if ~exist('max_iter')
    max_iter = 100;
end

if ~exist('assignment_flag')
    assignment_flag = 0;
end

pFlag = 1;
if ~exist('p')
    pFlag = 0;
else if p == 0;
        pFlag = 0;
    end
end

lambdaFlag = 1;
if ~exist('lambda')
    lambdaFlag = 0;
else if lambda == 0;
        lambdaFlag = 0;
    end
end

%--------------------------------------------------------------------------
%    Loop of Repetitions
%--------------------------------------------------------------------------

best = -1e100;
for itcount = 1:noruns

    fprintf('Iter. %g...  ',itcount);
    tic;

    %   initialize the parameters ------------------------------------------

    if pFlag == 0
        p = ones(1,k)/k;
    end

    if lambdaFlag == 0
        lambda = 500 ;
    end

    if mInitFlag == 0

        %
        %rand('seed',sum(100*clock));
        init = ceil(n*rand(k,1));
        mInit = x(init, :);
        %}
        %{
        v1init = ones(1,d)/sqrt(d);
        maskFarV1 = find(x*v1init' < 0.5);

        init = randselec(length(maskFarV1),k-1);
        mInit = [v1init; x(maskFarV1(init),:)];
        %}
    end

    m = mInit';
    mtemp = m;

    % ---------------------------------------------------------------------
    %   Main Loop ---------------------------------------------------------

    condition = 0;
    likelihood = [];
    temp = ones(1,n);

    if flagCorStop == 1

        mtcXitCor = [];

    end


    iter = 0;
    while condition == 0

        iter = iter + 1;

        % Estimation step   ---------------------------------------------------------

        dis = -(x*m) * lambda  - Cdln(lambda,d) ;
        mindis = min(dis,[],2);
        r = bsxfun(@times, p, exp(-bsxfun(@minus, dis, mindis)) );

        likelis = sum(r,2);
        r = bsxfun(@times, r, 1./likelis);
        p = sum(r,1)/n;

        % Maximization step -------------------------------------------------

        % Mean


        m = x'*r;
        m = bsxfun(@times, m, 1./sqrt(sum(m.^2,1)) );

        % Spread

        %{
          [lambdaNew fval exitflag]  = fzero(@(y) besseli(alpha+1,y)/besseli(alpha,y)-sum(sum(r.*(x*m)))/n,lambda);
          if exitflag == 1
              lambda = lambdaNew;
          end
        %}

        lambda = invAd(d,sum(sum(r.*(x*m)))/n);

        % Compute the likelihood

        likelihoodFinal = (1/n)*sum(log(likelis)-mindis);
        likelihood = [likelihood likelihoodFinal];



        % Checking the loop-break conditions ------------------------------
        if flagCorStop == 1

            % check parameters converge
            mtcXitCor = [mtcXitCor diag(m'*mtemp)];
            controlCondition = (sum(1-mtcXitCor(:,end) < epsilon) < k);
            mtemp = m;

        else

            if epsilonCheck == 0
                % check assignments converge
                [temp1 hardAssignments] = max(r');
                controlCondition = norm( hardAssignments - temp);
                temp = hardAssignments;

            else

                % check likelihood converge
                controlCondition = abs(likelihoodFinal - temp) > epsilon;
                temp = likelihoodFinal;

            end

        end

        % check assignments converge
        if(assignment_flag)
            [temp1 hardAssignments] = max(r');
            controlCondition = controlCondition + norm( hardAssignments - temp);
            temp = hardAssignments;
        end


        if (controlCondition <1)
            condition =1;
        end

        if (iter >= max_iter)
            condition = 1;
            warning(['Von mises did not converge after ' num2str(max_iter) ' iterations']);
        end
    end


    % Checking if this repetition gives the best solution -----------------

    [temp1 hardAssignments] = max(r');
    avcor = mean(sum(x.*m(:,hardAssignments)',2));

    if likelihoodFinal > best

        best = likelihoodFinal;
        clustered.lambda = lambda;
        clustered.likelihood = likelihood;
        clustered.r=r;
        clustered.clusters = hardAssignments';
        clustered.p = p;
        clustered.mtc = m';
        clustered.mtcInit = mInit;
        clustered.avcor = avcor;

        if flagCorStop == 1

            clustered.mtcXitCor = mtcXitCor;

        end

    end

    clustered.listlik = [clustered.listlik likelihoodFinal];
    clustered.listavcor = [clustered.listavcor avcor];
    timetemp = toc;
    fprintf('took %f second.\n',timetemp);

end

end

%-----------------------------------------------------------------------
function [out exitflag] = invAd(D,rbar)

outu = (D-1)*rbar/(1-rbar^2) + D/(D-1)*rbar;

[i] = besseli(D/2-1,outu);

if ((i == Inf)||(isnan(i)))
    out = outu - D/(D-1)*rbar/2;
    exitflag = Inf;
else
    [outNew fval exitflag]  = fzero(@(argum) Ad(argum,D)-rbar,outu);
    if exitflag == 1
        out = outNew;
    else
        out = outu - D/(D-1)*rbar/2;
    end
end
end
%-----------------------------------------------------------------------
function out = Ad(in,D)
% out = Ad(in,D)
out = besseli(D/2,in) ./ besseli(D/2-1,in);
end
%-----------------------------------------------------------------------
function out = Cdln(k,d)

% Computes the logarithm of the partition function of vonMises-Fisher as
% a function of kappa

sizek = size(k);
k = k(:);

out = (d/2-1).*log(k)-log(besseli((d/2-1)*ones(size(k)),k));

k0 = 500;
fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
nGrids = 1000;

maskof = find(k>k0);
nkof = length(maskof);

% The kappa values higher than the overflow

if nkof > 0

    kof = k(maskof);

    ofintv = (kof - k0)/nGrids;
    tempcnt = (1:nGrids) - 0.5;
    ks = k0 + bsxfun(@times, tempcnt, ofintv);
    adsum = sum( 1./((0.5*(d-1)./ks) + sqrt(1+(0.5*(d-1)./ks).^2)) ,2);

    % with correction:
    % 1./((0.5*(d-1)./x) + sqrt(1+(0.5*(d-1)./x).^2*(d/2/(d-1)^2)))

    out(maskof) =  fk0 - ofintv .* adsum;

end

out = reshape(out,sizek);
end

