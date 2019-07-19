function [flag, m2, SE, Va, Ve, Lnew] = Morphometricity_1RandEffect(y, X, K, alg, tol, MaxIter)
verbal = 0;

% This function implements morphometricity.
% It fits a linear mixed effects model using the restricted maximum likelihood (ReML) algorithm, 
% and produces the morphometricity estimate and its standard error. 
%
% The LME is: y = Xb + a + e, where X is the design (covariate) matrix of
% fixed effects, b is the fixed effect weight vector, a is a random effect
% drawn from a zero-mean Gaussian with a covariance equal to a scaled
% version of K (the ASM - see below), and e is i.i.d. zero mean Gaussian
% noise.

% Input -
% y is an Nsubj x 1 vector of phenotypes (trait values)
% X is an Nsubj x Ncov covariate matrix (that contains "nuisance variables
% such as age, sex, site dummy, etc)
% K is an Nsubj x Nsubj anatomical similarity matrix (ASM)
%  K has to be a symmetric, positive semi-definite matrix with its diagonal
%  elements averaging to 1.
%   If K is not non-negative definite, its negative eigenvalues will be set to zero,
%   and a warning will be printed to the command window
% alg is the algorithm for the ReML; default alg = 0
%   alg = 0 - use the average information
%   alg = 1 - use the expected information (Fisher scoring)
%   alg = 2 - use the observed information
% tol is the tolerance for the convergence of the ReML algorithm; default tol = 1e-4
% MaxIter is the maximum number of iterations for the ReML algorithm; default MaxIter = 100
%
% Output -
% flag indicates the convergence of the ReML algorithm
%   flag = 1 - the ReML algorithm has converged
%   flag = 0 - the ReML algorithm did not converged
% m2 is the morphometricity estimate
% SE is the standard error of the morphometricity estimate
% Va is the total anatomical/morphological variance
% Ve is the residual variance
% Lnew is the ReML likelihood when the algorithm is converged

%%%% Author: Tian Ge (minor modifications by Mert R. Sabuncu)
%%%% Contact: <tge1@nmr.mgh.harvard.edu> or <msabuncu@nmr.mgh.harvard.edu>

%% input check
if nargin < 3
    error('Not enough input arguments')
elseif nargin == 3
	alg = 0; tol = 1e-4; MaxIter = 100;
elseif nargin == 4
    tol = 1e-4; MaxIter = 100;
elseif nargin == 5
    MaxIter = 100;
elseif nargin > 6
    error('Too many input arguments')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit by Jingwei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(y,2)>1)
    error('''y'' should be a column vector.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit by Ruby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add bias term for X
X = [ones(size(y,1),1) X];

%check y in case there is nan
nan_index = find(isnan(y)==1);
y(nan_index) = [];
X(nan_index,:) = [];
K(nan_index,:) = [];
K(:,nan_index) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('In total %d subjects...\n', length(y));

[U, D] = eig(K);   % calculate the eigenvalues and eigenvectors of the GRM
if min(diag(D)) < 0   % check whether the GRM is non-negative definite
    disp('WARNING: the GRM is not non-negative definite! Set negative eigenvalues to zero')
    D(D<0) = 0;   % set negative eigenvalues to zero
    K = U*D/U;   % reconstruct the GRM
end
%% derived quantities
Nsubj = length(y);   % calculate the total number of subjects
Vp = var(y);   % calculate the phenotypic variance
%% initialization
Va = Vp/2; Ve = Vp/2;   % initialize the anatomical variance and residual variance
V = Va*K+Ve*eye(Nsubj);   % initialize the covariance matrix
P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;   % initialize the projection matrix
%% EM algorithm
if verbal
    disp('---------- EM algorithm ----------')
end
% use the expectation maximization (EM) algorithm as an initial update
Va = (Va^2*y'*P*K*P*y+trace(Va*eye(Nsubj)-Va^2*P*K))/Nsubj;   % update the anatomical variance
Ve = (Ve^2*y'*P*P*y+trace(Ve*eye(Nsubj)-Ve^2*P))/Nsubj;   % update the residual variance

% set negative estimates of the variance component parameters to Vp*1e-6
if Va < 0; Va = 10^(-6)*Vp; end 
if Ve < 0; Ve = 10^(-6)*Vp; end

V = Va*K+Ve*eye(Nsubj);   % update the covariance matrix
P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;   % update the projection matrix

E = eig(V); logdetV = sum(log(E));   % calculate the log determinant of the covariance matrix

Lold = inf; Lnew = -1/2*logdetV-1/2*log(det(X'/V*X))-1/2*y'*P*y;   % initialize the ReML likelihood
%% ReML
if verbal
disp('---------- ReML iterations ----------')
end
iter = 0;   % initialize the total number of iterations
while abs(Lnew-Lold)>=tol && iter<MaxIter   % criteria of termination
    
    % new iteration
    iter = iter+1; Lold = Lnew;
    if verbal
    disp(['---------- REML Iteration-', num2str(iter), ' ----------'])
    end
    % update the first-order derivative of the ReML likelihood
    Sg = -1/2*trace(P*K)+1/2*y'*P*K*P*y;   % score equation of the anatomical variance
    Se = -1/2*trace(P)+1/2*y'*P*P*y;   % score equation of the residual variance
    S = [Sg;Se];   % construct the score vector
    
    % update the information matrix
    if alg == 0
        Igg = 1/2*y'*P*K*P*K*P*y; Ige = 1/2*y'*P*K*P*P*y; Iee = 1/2*y'*P*P*P*y;   % average information
    elseif alg == 1
        Igg = 1/2*trace(P*K*P*K); Ige = 1/2*trace(P*K*P); Iee = 1/2*trace(P*P);   % expected information
    elseif alg == 2
        Igg = -1/2*trace(P*K*P*K)+y'*P*K*P*K*P*y; Ige = -1/2*trace(P*K*P)+y'*P*K*P*P*y; Iee = -1/2*trace(P*P)+y'*P*P*P*y;   % observed information
    end
    I = [Igg,Ige;Ige,Iee];   % construct the information matrix
    
    T = [Va;Ve]+I\S; Va = T(1); Ve = T(2);   % update the variance component parameters
    
    % set negative estimates of the variance component parameters to Vp*1e-6
    if Va < 0; Va = 10^(-6)*Vp; end 
    if Ve < 0; Ve = 10^(-6)*Vp; end
    
    V = Va*K+Ve*eye(Nsubj);   % update the covariance matrix
    P = (eye(Nsubj)-(V\X)/(X'/V*X)*X')/V;   % update the projection matrix
    
    if (isnan(sum(V(:))))
        flag = 0; m2 = NaN; SE = NaN; Va = NaN; Ve = NaN; Lnew = NaN; return;
    end
    E = eig(V); logdetV = sum(log(E+eps));   % calculate the log determinant of the covariance matrix
    
    Lnew = -1/2*logdetV-1/2*log(det(X'/V*X))-1/2*y'*P*y;   % update the ReML likelihood

end
%% morphometricity estimate and standard error
m2 = Va/(Va+Ve);   % morphometricity estimate

% update the information matrix at the final estimates
if alg == 0
    Igg = 1/2*y'*P*K*P*K*P*y; Ige = 1/2*y'*P*K*P*P*y; Iee = 1/2*y'*P*P*P*y;   % average information
elseif alg == 1
    Igg = 1/2*trace(P*K*P*K); Ige = 1/2*trace(P*K*P); Iee = 1/2*trace(P*P);   % expected information
elseif alg == 2
    Igg = -1/2*trace(P*K*P*K)+y'*P*K*P*K*P*y; Ige = -1/2*trace(P*K*P)+y'*P*K*P*P*y; Iee = -1/2*trace(P*P)+y'*P*P*P*y;   % observed information
end
I = [Igg,Ige;Ige,Iee];   % construct the score vector and the information matrix

invI = inv(I);   % inverse of the information matrix
SE = sqrt((m2/Va)^2*((1-m2)^2*invI(1,1)-2*(1-m2)*m2*invI(1,2)+m2^2*invI(2,2)));   % standard error estimate
%% diagnosis of convergence
if iter == MaxIter && abs(Lnew-Lold)>=tol
    flag = 0;
else
    flag = 1;
end
%%