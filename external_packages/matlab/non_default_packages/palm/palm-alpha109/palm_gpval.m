function pvals = palm_gpval(G,df1,df2,twotail)
% Convert a pivotal statistic computed with 'pivotal.m'
% (or simplifications) to a parametric p-value.
%
% Usage:
% pvals = palm_gpval(G,df1,df2)
%
% Inputs:
% G        : G or Z statistic.
% df1, df2 : Degrees of freedom (non infinite).
%            df1 must be a scalar.
%            For z, use df1 = 0.
%            For Chi2, use df1 = -1, and df2 as the df.
%
% Outputs:
% pvals    : Parametric p-values based on a
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
    
    % G or F, via conversion to Beta
    df2   = bsxfun(@times,ones(size(G)),df2);
    B     = (df1.*G./df2)./(1+df1.*G./df2);
    pvals = betainc(1-B,df2/2,df1/2);
    
elseif df1 == 1,
    
    % Student's t, Aspin-Welch v
    pvals = nan(size(G));
    Gsq = G.^2;
    df2 = bsxfun(@times,ones(size(G)),df2);
    ic  = df2 == 1;  % cauchy case
    in  = df2 > 1e7; % normal case
    ig  = ~(ic|in);  % general case
    if palm_isoctave,
        if any(ig),
            pvals(ig) = betainc(1./(1+Gsq(ig)./df2(ig)),df2(ig)/2,.5)/2;
        end
    else
        is  = df2 < Gsq;
        idx = ig(:) & is(:);
        if any(idx),
            pvals(idx) = betainc(1./(1+Gsq(idx)./df2(idx)),df2(idx)/2,.5)/2;
        end
        idx = ig(:) & ~is(:);
        if any(idx),
            pvals(idx) = betainc(1./(1+df2(idx)./Gsq(idx)),.5,df2(idx)/2,'upper')/2;
        end
    end
    ig  = G < 0 & ig;
    pvals(ig) = 1 - pvals(ig);
    if any(ic(:)),
        pvals(ic) = .5 + atan(-G(ic))/pi;
    end
    if any(in(:)),
        pvals(in) = palm_gpval(G(in),0);
    end
    pvals = reshape(pvals,size(G));
    
elseif df1 == 0,
    
    % Normal distribution
    pvals = erfc(G/sqrt(2))/2;
    
elseif df1 < 0,
    
    % Chi^2, via upper Gamma incomplete for precision and speed
    %df2   = bsxfun(@times,ones(size(G)),df2);
    pvals = palm_gammainc(G/2,df2/2,'upper');
    
end
