function gcdf = palm_gcdf(G,df1,df2)
% Convert a pivotal statistic computed with 'pivotal.m'
% (or simplifications) to a parametric p-value.
% The output is 1-p, i.e. the CDF.
% 
% Usage:
% cdf = palm_gcdf(G,df1,df2)
% 
% Inputs:
% G        : G or Z statistic.
% df1, df2 : Degrees of freedom (non infinite).
%            df1 must be a scalar
%            For z, use df1 = 0.
%            For Chi2, use df1 = -1, and df2 as the df.
% 
% Outputs:
% cdf      : Parametric cdf (1-p), based on a
%            t, F, z or Chi2 distribution.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Aug/2013
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

% Note that for speed, there's no argument checking,
% and some lines are repeated inside the conditions.

if df1 > 1,
    
    % G or F
    df2 = bsxfun(@times,ones(size(G)),df2);
    B = (df1.*G./df2)./(1+df1.*G./df2);
    gcdf = betainc(B,df1/2,df2/2);

elseif df1 == 1,
    
    % Student's t, Aspin's v
    df2 = bsxfun(@times,ones(size(G)),df2);
    ic = df2 == 1;
    in = df2 > 1e7;
    ig = ~(ic|in);
    gcdf = zeros(size(G));
    if any(ig(:)),
        gcdf(ig) = betainc(1./(1+G(ig).^2./df2(ig)),df2(ig)/2,.5)/2;
    end
    ig = G > 0 & ig;
    gcdf(ig) = 1 - gcdf(ig);
    if any(ic(:)),
        gcdf(ic) = .5 + atan(G(ic))/pi;
    end
    if any(in(:)),
        gcdf(ic) = erfc(-G(in)/sqrt(2))/2;
    end

elseif df1 == 0,
    
    % Normal distribution
    gcdf = erfc(-G/sqrt(2))/2;
    
elseif df1 < 0,
    
    % Chi^2, via lower Gamma incomplete for precision and speed
    %df2 = bsxfun(@times,ones(size(G)),df2);
    gcdf = palm_gammainc(G/2,df2/2,'lower');
    
end
