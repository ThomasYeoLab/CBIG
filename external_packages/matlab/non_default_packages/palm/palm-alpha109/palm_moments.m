function varargout = palm_moments(varargin)
% For a statistic that can be expressed as trace(A*W), for
% a sample size of n observations, this function returns the
% expected first three moments of the permutation distribution,
% without actually computing any permutation.
%
% [mu,sigsq,gamm1,gamm2] = palm_moments(G)
% [mu,sigsq,gamm1]       = palm_moments(A,W,n)
%
% Inputs:
% - G : A PxV array of observations of the G random variable.
%       The moments are unbiased and run along the 1st dimension.
%       Typical case is P = number of permutations and V = 
%       number of tests, e.g., voxels.
% - A : A square matrix for a multivariate proper statistic,
%       or a vector of A values for various univariate tests.
% - W : A square matrix for a multivariate proper statistic,
%       or a vector of W values for various univariate tests.
% - n : Sample size on which the statistic is based.
%
% Outputs:
% - mu    : Sample mean.
% - sigsq : Sample variance (unbiased).
% - gamm1 : Sample skewness (unbiased).
% - gamm2 : Sample kurtosis (unbiased).
%
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% For the estimators using trace(AW), the references are:
% * Kazi-Aoual F, Hitier S, Sabatier R, Lebreton J-D. Refined
%   approximations to permutation tests for multivariate
%   inference. Comput Stat Data Anal. 1995;20(94):643-656.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;4:497-511.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Mar/2015
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin == 1,
    
    % For a set of values of the random variable G, return the
    % first 4 moments.
    
    % Mean
    G     = varargin{1};
    n     = size(G,1);
    mu    = sum(G,1)/n;
    
    % Variance
    G0    = bsxfun(@minus,G,mu);
    ssq   = sum(G0.^2,1);
    sigsq = (ssq/(n-1));
    
    % Skewness
    s2    = ssq/n; % biased variance
    m3    = sum(G0.^3,1)/n;
    gamm1 = m3./s2.^1.5;
    gamm1 = gamm1 * sqrt((n-1)/n)*n/(n-2); % unbiased skewness

    % Kurtosis (normal dist = 3)
    if nargout == 4,
        m4    = sum(G0.^4,1)/n;
        gamm2 = (m4./s2.^2);
        gamm2 = ((n+1)* gamm2 -3*(n-1))*(n-1)/((n-2)*(n-3))+3; % unbiased kurtosis
    else
        gamm2 = [];
    end
    
elseif nargin == 3,
    
    % Compute the first three moments of the permutation distribution of
    % the statistic G = trace(AW), for n subjects, using the method in
    % Kazi-Aoual et al (1995). The variable names follow roughly the same
    % as in the paper.
    
    % Take inputs
    A = varargin{1};
    W = varargin{2};
    n = varargin{3};
    
    % If A and W are truly multivariate (i.e., square matrices), do as in
    % the original paper. Otherwise, make simplifications as these are all
    % scalars.
    if size(A,1) == size(A,2),
        
        % Some auxiliary variables for ET2:
        T    = trace(A);
        T_2  = trace(A^2);
        S_2  = sum(diag(A).^2);
        Ts   = trace(W);
        Ts_2 = trace(W^2);
        Ss_2 = sum(diag(W).^2);
        
        % Some auxiliary variables for ET3:
        T_3  = trace(A^3);
        S_3  = sum(diag(A).^3);
        U    = sum(A(:).^3);
        R    = diag(A)'*diag(A^2);
        B    = diag(A)'*A*diag(A);
        Ts_3 = trace(W^3);
        Ss_3 = sum(diag(W).^3);
        Us   = sum(W(:).^3);
        Rs   = diag(W)'*diag(W^2);
        Bs   = diag(W)'*W*diag(W);
        
    else
        
        % Some auxiliary variables for ET2:
        T    = A;
        T_2  = A.^2;
        S_2  = T_2;
        Ts   = W;
        Ts_2 = W.^2;
        Ss_2 = Ts_2;
        
        % Some auxiliary variables for ET3:
        T_3  = A.^3;
        S_3  = T_3;
        U    = T_3;
        R    = T_3;
        B    = T_3;
        Ts_3 = W.^3;
        Ss_3 = Ts_3;
        Us   = Ts_3;
        Rs   = Ts_3;
        Bs   = Ts_3;
    end
    
    % E(T):
    mu = T.*Ts/(n-1);
    
    % V(T):
    sigsq = 2*((n-1)*T_2-T.^2).*((n-1)*Ts_2-Ts.^2) / (n-1)^2/(n+1)/(n-2) ...
        + (n*(n+1)*S_2-(n-1)*(T.^2+2*T_2)) .* (n*(n+1)*Ss_2-(n-1)*(Ts.^2+2*Ts_2)) ...
        / (n+1)/n/(n-1)/(n-2)/(n-3);
    
    % E(T^3):
    ET3 = ...
        n^2*(n+1)*(n^2+15*n-4)*S_3.*Ss_3 ...
        + 4*(n^4-8*n^3+19*n^2-4*n-16)*U.*Us ...
        + 24*(n^2-n-4)*(U.*Bs+B.*Us) ...
        + 6*(n^4-8*n^3+21*n^2-6*n-24)*B.*Bs ...
        + 12*(n^4-n^3-8*n^2+36*n-48)*R.*Rs ...
        + 12*(n^3-2*n^2+9*n-12)*(T.*S_2.*Rs + R.*Ts.*Ss_2) ...
        + 3*(n^4-4*n^3-2*n^2+9*n-12)*T.*Ts.*S_2.*Ss_2 ...
        + 24*( (n^3-3*n^2-2*n+8)*(R.*Us+U.*Rs) ...
        + (n^3-2*n^2-3*n+12)*(R.*Bs+B.*Rs) ) ...
        + 12*(n^2-n+4)*(T.*S_2.*Us+U.*Ts.*Ss_2) ...
        + 6*(2*n^3-7*n^2-3*n+12)*(T.*S_2.*Bs+B.*Ts.*Ss_2) ...
        - 2*n*(n-1)*(n^2-n+4)*( (2*U+3*B).*Ss_3+(2*Us+3*Bs).*S_3 ) ...
        - 3*n*(n-1)^2*(n+4)*( (T.*S_2+4*R).*Ss_3+(Ts.*Ss_2+4*Rs).*S_3 ) ...
        + 2*n*(n-1)*(n-2)*( (T.^3+6*T.*T_2+8*T_3).*Ss_3 ...
        + (Ts.^3+6*Ts.*Ts_2+8*Ts_3).*S_3 ) ...
        + T.^3.*((n^3-9*n^2+23*n-14)*Ts.^3+6*(n-4).*Ts.*Ts_2+8*Ts_3) ...
        + 6*T.*T_2.*((n-4)*Ts.^3+(n^3-9*n^2+24*n-14)*Ts.*Ts_2+4*(n-3)*Ts_3) ...
        + 8*T_3.*(Ts.^3+3*(n-3).*Ts.*Ts_2+(n^3-9*n^2+26*n-22)*Ts_3) ...
        - 16*(T.^3.*Us+U.*Ts.^3)-6*(T.*T_2.*Us+U.*Ts.*Ts_2)*(2*n^2-10*n+16) ...
        - 8*(T_3.*Us+U.*Ts_3)*(3*n^2-15*n+16)-(T.^3.*Bs+B.*Ts.^3) ...
        * (6*n^2-30*n+24)-6*(T.*T_2.*Bs+B.*Ts.*Ts_2)*(4*n^2-20*n+24) ...
        - 8*(T_3.*Bs + B.*Ts_3)*(3*n^2-15*n+24) ...
        - (n-2)*( 24*(T.^3.*Rs+R.*Ts.^3)+6*(T.*T_2.*Rs+R.*Ts.*Ts_2)*(2*n^2-10*n+24) ...
        + 8*(T_3.*Rs+R.*Ts_3)*(3*n^2-15*n+24)+(3*n^2-15*n+6) ...
        .* (T.^3.*Ts.*Ss_2+T.*S_2.*Ts.^3) ...
        + 6*(T.*T_2.*Ts.*Ss_2+T.*S_2.*Ts.*Ts_2)*(n^2-5*n+6) ...
        + 48*(T_3.*Ts.*Ss_2+T.*S_2.*Ts_3) );
    ET3 = ET3/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5);
    
    % The coefficient "3" below is missing from Kazi-Aoual et al (1995), but it
    % is shown in the Supplementary Information of Minas and Montana (2014).
    gamm1 = (ET3 - 3*mu.*sigsq - mu.^3)./sigsq.^1.5;
    gamm2 = [];
else
    error('Incorrect number of arguments');
end

% Return results
varargout{1} = mu;
varargout{2} = sigsq;
varargout{3} = gamm1;
varargout{4} = gamm2;
