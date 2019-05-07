function Z = palm_gtoz(G,df1,df2)
% Convert a G-statistic (or any of its particular cases)
% to a z-statistic (normally distributed).
%
% Usage:
% Z = palm_gtoz(G,df1,df2)
%
% Inputs:
% G        : G statistic.
% df1, df2 : Degrees of freedom (non-infinite).
% 
% Outputs:
% Z        : Z-score
%
% If df2 = NaN and df1 = 1, G is treated as Pearson's r.
% If df2 = NaN and df1 > 1, G is treated as R^2.
% If df2 = NaN and df1 = 0, G is treated as z already.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Jan/2014
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

% Note that for speed, there's no argument checking.

% If df2 is NaN, this is r, R^2, or z already
if isnan(df2(1)),
    
    if df1 == 0,
        
        % If df1 is zero, this is already a z-stat (this is here more for
        % compatibility).
        Z = G;
        
    elseif df1 == 1,
        
        % If rank(C) = 1, i.e., df1 = 1, this is r, so
        % do a Fisher's r-to-z stransformation
        Z = atanh(G);
        
    elseif df1 > 1,
        
        % If rank(C) > 1, i.e., df1 > 1, this is R^2, so
        % use a probit transformation.
        Z = -erfcinv(2*G)*sqrt(2); %Z = norminv(G);
        
    end

else
    siz = size(G);
    Z   = zeros(siz);
    df2 = bsxfun(@times,ones(siz),df2);
    if df1 == 1,
        
        % Deal with precision issues working on each
        % tail separately
        idx = G > 0;
        %Z( idx) = -erfinv(2*palm_gcdf(-G( idx),1,df2( idx))-1)*sqrt(2);
        %Z(~idx) =  erfinv(2*palm_gcdf( G(~idx),1,df2(~idx))-1)*sqrt(2);
        Z( idx) =  erfcinv(2*palm_gcdf(-G( idx),1,df2( idx)))*sqrt(2);
        Z(~idx) = -erfcinv(2*palm_gcdf( G(~idx),1,df2(~idx)))*sqrt(2);
        
    elseif df1 == 0,
        
        % If df1 is zero, this is already a z-stat (this is here more for
        % compatibility).
        Z = G;
        
    else
        
        % G-vals above the upper half are treated as
        % "upper tail"; otherwise, "lower tail".
        thr = (1./betainv(.5,df2/2,df1/2)-1).*df2/df1;
        idx = G > thr;
        
        % Convert G-distributed variables to Beta-distributed
        % variables with parameters a=df1/2 and b=df2/2
        B = (df1.*G./df2)./(1+df1.*G./df2);
        a = df1/2;
        b = df2/2;
        
        % Convert to Z through a Beta incomplete function
        %Z( idx) = -erfinv(2*betainc(1-B( idx),b( idx),a)-1)*sqrt(2);
        %Z(~idx) =  erfinv(2*betainc(  B(~idx),a,b(~idx))-1)*sqrt(2);
        Z( idx) =  erfcinv(2*betainc(1-B( idx),b( idx),a))*sqrt(2);
        Z(~idx) = -erfcinv(2*betainc(  B(~idx),a,b(~idx)))*sqrt(2);
        
    end
end