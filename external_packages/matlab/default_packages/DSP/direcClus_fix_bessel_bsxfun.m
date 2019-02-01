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
% !!IMPORTANT!!
%       Please note that in this algorithm the mean directions of clusters
%       will be initialized by random seeds. Therefore we will use a rand
%       function to set the random seeds. Please set your own random number
%       generator before you run this script, so that you can obtain the
%       same results if you re-run this script. To do that, you can use
%       rng(seed) to seeds the random number generator (see 
%       https://www.mathworks.com/help/matlab/ref/rng.html for more detail). 
%       
%       **Parallel random initializations**
%       In our current code, if you set noruns = 1000, it will sequentially
%       run 1000 times and pick the best result. This should be fine if 
%       each initialization is very fast, however, if you find each 
%       initialization is quite slow, you may consider run multiple random 
%       initializations in parallel. For example, you can split 1000 tries 
%       as 100 jobs, where each job will run a 10 tries clustering. To do 
%       this, you just need to give different random seed for each job, 
%       i,e, you can do rng(1) for job1,  and rng(2) for job2. After you 
%       obtain the results of 100 jobs, then you just need to pick the one 
%       with the highest clustered.likelihood(end) as your final output.
%
% INPUT: 
%
%       - x
%         NxD matrix. The input data. N is the number of vertices, D is 
%         the number of features.
%
%       - k
%         A scalar. The number of clusters you want to perform the 
%         clustering. For example, in Yeo2011, we do 17 or 7.
%
%       - d
%         A scalar. The dimensionality of features - 1, i.e. d = size(x,2)-1.
%
%       - noruns
%         A scalar. The number of random initializations. When we do 
%         clustering, the algorithm starts with random seeds, if you set 
%         noruns = 1000, it will run 1000 different random seeds and pick 
%         the most optimal one as the final output.
%
%       - lambda
%         A scalar. This is the lowest variance that allows for each 
%         cluster. We will use it to for initialization. Noted that lower 
%         lambda indicates higher variance in vmf, therefore lambda should 
%         be small. If lambda = 0, the script will automatically estimate it.
%
%       - p (set p = 0 if no special reason)
%         A 1xk vector. The initialization parameter of network
%         probability for each vertex. We recommend the user to set p = 0.
%         Then the default p = ones(1,k)/k, which indicates equal
%         probability of assigning a vertex into a cluster. k is the number
%         of clusters. If p is given, the code will use the given p for
%         initialization.
%
%       - mInit (set mInit = 0 if no special reason)
%         A kxD matrix. The initial seeds of networks. We recommend the
%         user to set mInit = 0, then the clustering will be performed with
%         k random seeds, the dimensionality of each seed is D. 
%
%
%       - epsilon (set epsilon = 1e-4 if no special reason)
%         A scalar. A small value for checking convergence.
%
%       - flagCorStop (set flagCorStop = 1 if no special reason)
%         1 or 0. The stoppig criteria of convergence for EM iterations. If
%         flagCorStop = 1, will check whether the mean direction of vmf
%         algorithm can converge, i.e. it uses the correlation between the 
%         updated mean directions and the old ones in each iteration as the
%         condition for breaking the loop, and does so by having such 
%         correlations for all the k mean vectors to be at least
%         1-epsilon. If flagCorStop = 0, will check whether the likelihood
%         cost can converge, i.e. it uses the changes between the updated
%         cost and old cost in each iteration as the condition for breaking
%         the loop.
%
%       - max_iter (set max_iter = 100 if no speacial reason)
%         A scalar. The maximum number of EM iterations.
%
%       - assignment_flag (set assignment_flag = 1 if no special reason)
%         1 or 0. If assignment_flag = 1, the algorithm will check whether
%         the cluster assignment can converge.
%
% OUTPUT:
%
%       - clustered 
%         A structure, there are some fields that are important: 
%
%       1. clustered.clusters: 
%          It should be a Nx1 vector contains the final clustering results. 
%          For example, if N=4, num_clusters = 2, and 
%          clustered.clusters = [1;2;2;1]. This means vertex 1 and 4 are 
%          clustered together and vertex 2 and 3 are clustered together. 
%          Note that clustered.clusters = [2;1;1;2] is the same result as 
%          clustered.clusters= [1;2;2;1].
%
%       2. clustered.likelihood: 
%          It should be a 1 x num_tries vector contains the cost for each 
%          initialization. The last value, i.e, clustered.likelihood(end) 
%          is the final cost. (The higher the cost the better the performance).
%
%       3. clustered.lambda: 
%          A scalar. The estimated concentration parameter (higher lambda 
%          indicates lower variance)
%
%       4. clustered.mtc: 
%          A num_clusters x D matrix. The estimated mean direction for each
%          cluster, each row of clustered.mtc represent the mean direction 
%          for a cluster. You can interpret it as the average of all the 
%          vertices within each cluster. 
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
x = bsxfun(@minus, x, mean(x,2));

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
else
    if lambda == 0;
        lambdaFlag = 0;
    else
        inilambda = lambda;
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
        ini_flag = 0;
        start_value = 50;
        ini_loop = 0;
        while (ini_flag == 0)
            ini_loop = ini_loop + 1;
            [inilambda,~,exitflag] = fzero(@(inputx) sign(inputx)*abs(besseli(alpha,inputx))-1e+10, start_value,optimset('Display','off'));
            if (exitflag == 1)
                lambda = inilambda;
                ini_flag = 1;
            else
                start_value = start_value + 50;
%                 fprintf('try another starting value...%d\n',start_value);
            end
            if(ini_loop == 1000)
                error('Can not initialize lambda automatically, please manually set it');
            end
        end
        fprintf('initialize lambda with %.2f\n',lambda);
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

        dis = -(x*m) * lambda  - Cdln(lambda,d,inilambda) ;
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
function out = Cdln(k,d,k0)

% Computes the logarithm of the partition function of vonMises-Fisher as
% a function of kappa

sizek = size(k);
k = k(:);

out = (d/2-1).*log(k)-log(besseli((d/2-1)*ones(size(k)),k));

% find k0:

fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
nGrids = 1000;

if isinf(fk0)% more general function works at least up to d=350k
    k0=0.331*d;
    fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
end
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

