function [F,Fdf1,Fdf2] = palm_qtof(Q,df1,df2,p,Qname)
% Convert a multivariate MANOVA/MANCOVA scalar
% statistic (Q) to an F-score. The conversion is only
% approximate for most methods and cases.
% 
% Usage:
% [F,Fdf1,Fdf2] = palm_qtof(Q,df1,df2,p,Qname)
% 
% Inputs:
% - Q     : Q-statistic.
% - df1   : Degrees of freedom of the model. This is typically
%           the rank of the contrast, rank(C).
% - df2   : Degrees of freedom of the error. This is typically
%           N-rank(M).
% - p     : Number of dependent variables in the model. This is
%           typically rank(Y) or rank(D).
% - Qname : Name of the stastistic. It can be one of:
%           'Wilks', 'Pillai', 'Hotelling' or 'Roy_ii'.
% 
% Outputs:
% - F     : F-statistic (often approximate).
% - Fdf1  : Degrees of freedom of the hypothesis.
% - Fdf2  : Degrees of freedom of the error.
% 
% _____________________________________
% Anderson Winkler
% FMRIB / University of Oxford
% Mar/2014
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

switch lower(Qname),
    case 'wilks',
        r = df2-(p-df1+1)/2;
        u = (p*df1-2)/4;
        cden = (p^2+df1^2-5);
        if cden > 0,
            t = sqrt((p^2*df1^2-4)/cden);
        else
            t = 1;
        end
        F = (r*t-2*u)*(1-Q.^(1/t))./(Q.^(1/t)*p*df1);
        Fdf1 = p*df1;
        Fdf2 = r*t-2*u;

    case 'pillai',
        m = (abs(p-df1)-1)/2;
        n = (df2-p-1)/2;
        s = min(p,df1);
        F = (2*n+s+1)/(2*m+s+1)*(Q./(s-Q));
        Fdf1 = s*(2*m+s+1);
        Fdf2 = s*(2*n+s+1);
        
    case 'hotelling',
        m = (abs(p-df1)-1)/2;
        n = (df2-p-1)/2;
        s = min(p,df1);
        if n > 0,
            b = (p+2*n)*(df1+2*n)/(2*(2*n+1)*(n-1));
            c = (2+(p*df1+2)/(b-1))/(2*n);
            Fdf1 = p*df1;
            Fdf2 = 4+(p*df1+2)/(b-1);
            F = (Q/c)*Fdf2/Fdf1;
        else
            Fdf1 = s*(2*m+s+1);
            Fdf2 = 2*(s*n+1);
            F = (Q/s)*Fdf2/Fdf1;
        end
        
    case 'roy_ii',
        Fdf1 = max(p,df1);
        Fdf2 = df2-Fdf1+df1;
        F = Q*Fdf2/Fdf1;
        
    otherwise
        error('Statistic %s unknown.',Qname);
end