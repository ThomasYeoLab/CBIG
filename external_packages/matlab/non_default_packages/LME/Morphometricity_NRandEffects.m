function [flag, m2_tot, se_tot, M2, SE, Vc, Lnew] = Morphometricity_NRandEffects(y, X, K, alg, tol, max_iter, verbal)
% This function implements morphometricity for multiple kernels.
% It fits a linear mixed effects model using the restricted maximum likelihood (ReML) algorithm, 
% and produces morphometricity estimates and their standard errors. 
%
% Input -
% y is an Nsubj x 1 vector of phenotypes
% X is an Nsubj x Ncov covariate matrix
%   Note that intercept will NOT be automatically added to the model and need to be
%   explicitly included in the covariate matrix if needed
% K is an Nsubj x Nsubj x Nker matrix, comprising Nker similarity matrices
%   Note that an identity matrix will be automatically added to model the covariance of residuals (independent environment)
% alg is the algorithm for the ReML; default alg = 0
%   alg = 0 - use the average information
%   alg = 1 - use the expected information (Fisher scoring)
%   alg = 2 - use the observed information
% tol is the tolerance for the convergence of the ReML algorithm; default tol = 1e-4
% max_iter is the maximum number of iterations for the ReML algorithm; default max_iter = 100
% verbal = true - print number of iterations to the Command Window; default verbal = false
%
% Output -
% flag indicates the convergence of the ReML algorithm
%   flag = 1 - the ReML algorithm has converged
%   flag = 0 - the ReML algorithm did not converged
% m2_tot - total fraction of phenotypic variance explained by all the input kernels
% se_tot - standard error estimate of m2_tot
% M2 is an (Nker+1) x 1 vector of morphometricity estimates (proportion of phenotypic variance explained by each kernel)
% SE is an (Nker+1) x 1 vector of standard errors for the morphometricity estimates M2
% Vc is an (Nker+1) x 1 vector of the estimates of variance component parameters
% Lnew is the ReML likelihood when the algorithm is converged

%% input check
if nargin < 3
    error('Not enough input arguments')
elseif nargin == 3
	alg = 0; tol = 1e-4; max_iter = 100; verbal = false;
elseif nargin == 4
    tol = 1e-4; max_iter = 100; verbal = false;
elseif nargin == 5
    max_iter = 100; verbal = false;
elseif nargin == 6
    verbal = false;
elseif nargin > 7
    error('Too many input arguments')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit by Jingwei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(y, 2)>1)
    error('''y'' should be a column vector.');
end

% Add bias term for X
X = [ones(size(y,1),1) X];

%check y in case there is nan
nan_index = find(isnan(y)==1);
y(nan_index) = [];
X(nan_index,:) = [];
K(nan_index,:,:) = [];
K(:,nan_index,:) = [];



%% derived quantities
Nsubj = length(y);   % calculate the total number of subjects
Vp = var(y);   % calculate the phenotypic variance
K(:,:,end+1) = eye(Nsubj);   % add an identity matrix to model the covariance of residuals
Nk = size(K,3);   % total number of variance components
%% initialization
Vc = Vp/Nk*ones(Nk,1);   % initialize variance component parameters

% initialize the covariance matrix
V = zeros(Nsubj,Nsubj);
for i = 1:Nk
    V = V+Vc(i)*K(:,:,i);
end

P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;   % initialize the projection matrix
%% EM algorithm
if verbal
    disp('---------- EM algorithm ----------')
end
% use the expectation maximization (EM) algorithm as an initial update
for i = 1:Nk
    Vc(i) = (Vc(i)^2*y'*P*K(:,:,i)*P*y+trace(Vc(i)*eye(Nsubj)-Vc(i)^2*P*K(:,:,i)))/Nsubj;
end
Vc(Vc<0) = 10^(-6)*Vp;   % set negative estimates of the variance component parameters to Vp*1e-6

% update the covariance matrix
V = zeros(Nsubj,Nsubj);
for i = 1:Nk
    V = V+Vc(i)*K(:,:,i);
end

P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;   % update the projection matrix

E = eig(V); logdetV = sum(log(E));   % calculate the log determinant of the covariance matrix
Lold = inf; Lnew = -1/2*logdetV-1/2*log(det(X'/V*X))-1/2*y'*P*y;   % initialize the ReML likelihood
%% ReML
if verbal
    disp('---------- ReML iterations ----------')
end

iter = 0;   % initialize the total number of iterations
while abs(Lnew-Lold)>=tol && iter<max_iter   % criteria of termination
    
    % new iteration
    iter = iter+1; Lold = Lnew;
    if verbal
        disp(['---------- Iteration-', num2str(iter), ' ----------'])
    end
    % update the first-order derivative of the ReML likelihood
    S = zeros(Nk,1);
    for i = 1:Nk
        S(i) = -1/2*trace(P*K(:,:,i))+1/2*y'*P*K(:,:,i)*P*y;
    end
    
    % update the information matrix
    I = zeros(Nk,Nk);
    if alg == 0
        for i = 1:Nk
            for j = 1:Nk
                I(i,j) = 1/2*y'*P*K(:,:,i)*P*K(:,:,j)*P*y;   % average information
            end
        end
    elseif alg == 1
        for i = 1:Nk
            for j = 1:Nk
                I(i,j) = 1/2*trace(P*K(:,:,i)*P*K(:,:,j));   % expected information
            end
        end
    elseif alg == 2
        for i = 1:Nk
            for j = 1:Nk
                I(i,j) = -1/2*trace(P*K(:,:,i)*P*K(:,:,j))+y'*P*K(:,:,i)*P*K(:,:,j)*P*y;   % observed information
            end
        end
    end
    
    Vc = Vc+I\S;   % update the variance component parameters
    Vc(Vc<0) = 10^(-6)*Vp;   % set negative estimates of the variance component parameters to Vp*1e-6

    % update the covariance structure
    V = zeros(Nsubj,Nsubj);
    for i = 1:Nk
        V = V+Vc(i)*K(:,:,i);
    end
    
    P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;   % update the projection matrix
    
    E = eig(V); logdetV = sum(log(E));   % calculate the log determinant of the covariance matrix    
    Lnew = -1/2*logdetV-1/2*log(det(X'/V*X))-1/2*y'*P*y;   % update the ReML likelihood
end
%% diagnosis of convergence
if iter == max_iter && abs(Lnew-Lold)>=tol
    flag = 0;
else
    flag = 1;
end
%% morphometricity estimates and standard errors
m2_tot = sum(Vc(1:end-1))/sum(Vc);   % total fraction of phenotypic variance explained
M2 = Vc/sum(Vc);   % morphometricity estimates

% update the information matrix evaluated at the final estimates
I = zeros(Nk,Nk);
if alg == 0
    for i = 1:Nk
        for j = 1:Nk
            I(i,j) = 1/2*y'*P*K(:,:,i)*P*K(:,:,j)*P*y;   % average information
        end
    end
elseif alg == 1
    for i = 1:Nk
        for j = 1:Nk
            I(i,j) = 1/2*trace(P*K(:,:,i)*P*K(:,:,j));   % expected information
        end
    end
elseif alg == 2
    for i = 1:Nk
        for j = 1:Nk
            I(i,j) = -1/2*trace(P*K(:,:,i)*P*K(:,:,j))+y'*P*K(:,:,i)*P*K(:,:,j)*P*y;   % observed information
        end
    end
end

% standard error estimates
T = zeros(Nk,1);
for i = 1:Nk
    if i == Nk
        T(i) = -m2_tot/sum(Vc);
    else
        T(i) = M2(end)/sum(Vc);
    end
end
se_tot = sqrt(T'/I*T);

SE = zeros(Nk,1);
for i = 1:Nk
    T = zeros(Nk,1);
    for j = 1:Nk
        if j == i
            T(j) = M2(i)*(1-M2(i))/Vc(i);
        else
            T(j) = -M2(i)^2/Vc(i);
        end
    end
    SE(i) = sqrt(T'/I*T);
end
%%