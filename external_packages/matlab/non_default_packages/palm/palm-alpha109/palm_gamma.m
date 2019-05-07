function pvals = palm_gamma(G,mu,sigsq,gamm1,rev,prepl)
% Return the p-values for a Gamma distribution, parameterised by
% its first three moments.
%
% pvals = palm_gamma(G,mu,s2,gamm1,rev)
% 
% Inputs:
% - G     : Statistics for which p-values are to be computed.
% - mu    : Distribution mean.
% - sigsq : Distribution standard deviation.
% - gamm1 : Distribution skewness.
% - rev   : Use if lower values of the statistic are evidence in
%           favour of the alternative.
% - prepl : Replacement for what otherwise would be zero p-values
%           in case of poor fits (e.g., statistic falls into the
%           part of the distribution that has pdf=0. In these cases
%           the p-value can be 1 or 1/(#perm) depending on which
%           tail and the sign of the skewness.
%
% Outputs:
% - pvals : p-values.
% 
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% Other references:
% * Mielke PW, Berry KJ, Brier GW. Application of Multi-Response
%   Permutation Procedures for Examining Seasonal Changes in
%   Monthly Mean Sea-Level Pressure Patterns. Mon Weather Rev.
%   1981;109(1):120-126.
% * Minas C, Montana G. Distance-based analysis of variance:
%   Approximate inference. Stat Anal Data Min. 2014;7(6):450-470.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2015
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

% Note that there are no argument checking for speed, but
% sizes of all inputs need to be the same, or the moments need to
% be all scalars.

if gamm1 == 0,
    
    % If not skewed, use a normal approximation.
    G     = (G - mu)./sigsq.^.5;
    pvals = erfc(G/sqrt(2))/2;
    
else
    
    % Standardise G, so that all becomes a function of the skewness.
    G    = (G - mu)./sigsq.^.5;
    
    % Gamma distribution parameters (Minas & Montana, 2014).
    kpar = 4/gamm1.^2;
    tpar = gamm1/2;
    cpar = -2/gamm1;
     
    % Actual p-value. If there are negatives here, the probability can
    % have an imaginary part, which is dealt with later.
    if rev,
        if gamm1 > 0,
            tail = 'lower';
        else
            tail = 'upper';
        end
    else
        if gamm1 > 0,
            tail = 'upper';
        else
            tail = 'lower';
        end
    end
    pvals = palm_gammainc((G-cpar)./tpar,kpar,tail);
    
    % Deal with imaginary parts.
    if ~ isreal(pvals),
        iidx = imag(pvals) ~= 0;
        if rev,
            if gamm1 > 0,
                pvals(iidx) = prepl;
            else
                pvals(iidx) = 1;
            end
        else
            if gamm1 > 0,
                pvals(iidx) = 1;
            else
                pvals(iidx) = prepl;
            end
        end
    end
end
