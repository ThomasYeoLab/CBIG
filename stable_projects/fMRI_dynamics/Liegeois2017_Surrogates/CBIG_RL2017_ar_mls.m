function [Y,B,Z,E] = CBIG_RL2017_ar_mls(TS,p)

% [Y,B,Z,E] = CBIG_RL2017_ar_mls(TS,p)
% This function identifies an autoregressive (AR) model from multivariate time series using the multivariate least square approach

% Input

% 'TS'            original time course encoded in a matrix of size (T=number of time points,k=number of variables/ROIs)
% 'p'             order of the AR model to be identified

% Given these inputs, the AR model of order 'p' reads:

% TS(t,:)' = w + \sum_i=1^p A_i * TS(t-i,:)' + e(t),     Eq. (1)

% where: w    is an intercept vector of size (T,1). It is supposed to be
%             zero if timecourses are centered.
%        A_i  are matrices of size k*k linking T(t,:) and T(t-i,:).
%        e(t) is an error vector of size (T,1) following a centered multivariate gaussian distribution

% Eq. (1) can be written for all t \in [p+1...T]. Concatenating these
% (T-p) equations yields the following matrix form:

% Y = B*Z+E,    Eq. (2)

% where: Y = [TS(p+1,:)' TS(p+2,:)' ... TS(T,:)'] is a matrix of size (k,T-p)
%
%        B = [w A_1 ... A_p] is a matrix of size (k,k*p+1) that
%        gathers unknown parameters of the model
%
%            |    1             1           ...        1     |
%            | TS(p,:)'     TS(p+1,:)'      ...    TS(T-1,:)'|
%            |TS(p-1,:)'     TS(p,:)'       ...    TS(T-2,:)'|
%        Z = |    .             .          .           .     |
%            |    .             .            .         .     |
%            |    .             .              .       .     |
%            | TS(1,:)'      TS(2,:)'       ...    TS(T-p,:)'|
%
%        is a matrix of size (k*p+1,T-p) that is directly built from the input TS.
%
%        E = [e(p+1) e(p+2) ... e(T)] is a matrix of size (k,T-p)
%        containing the residuals of the multivariate AR model.


% Output
%
% 'Y'             Matrix variables directly built from TS (see Eq. (2))
% 'B'             Matrix containing AR model parameters
% 'Z'             Matrix variables directly built from TS (see Eq. (2))
% 'E'             Residuals of the AR model

% Note that the most important output is B which contains the parameters of
% the AR model. E is also of interest when autoregressive randomization is
% to be performed. Y and Z are provided for information and testing
% purposes.

% Written by Raphael Liegeois and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Contact: Raphael Liegeois (Raphael.Liegeois@gmail.com)



%% Identified parameters
T   = size(TS,1);
k   = size(TS,2);


% Rewriting the problem as Y=BZ+E
Y   = zeros(k,T-p);
Z   = zeros(k*p+1,T-p);

% Filling up Y from the data TS

Y = TS';
Y(:, 1:p) = [];

% Filling up Z

% First row (intercept)

Z(1,:)=ones(1,T-p);

% Other rows

for j = 1:p
    Z((j-1)*k+2 : j*k+1, :) = TS(p-j+1:T-j, :)';
end

% Identifying AR parameters

B   = (Y*Z')/(Z*Z');

% Computing residuals

E   = Y-B*Z;




