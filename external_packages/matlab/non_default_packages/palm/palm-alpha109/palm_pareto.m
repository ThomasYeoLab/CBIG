function [P,apar,kpar,upar] = palm_pareto(G,Gdist,rev,Pthr,G1out)
% Compute the p-values for a set of statistics G, taking
% as reference a set of observed values for G, from which
% the empirical cumulative distribution function (cdf) is
% generated. If the p-values are below Pthr, these are
% refined further using a tail approximation from the
% Generalised Pareto Distribution (GPD).
%
% Usage:
% P = palm_pareto(G,Gdist,rev,Pthr)
%
% Inputs:
% G      : Vector of Nx1 statistics to be converted to p-values
% Gdist  : A Mx1 vector of observed values for the same statistic
%          from which the empirical cdf is build and p-values
%          obtained. It doesn't have to be sorted.
% rev    : If true, indicates that the smallest values in G and
%          Gvals, rather than the largest, are the most significant.
% Pthr   : P-values below this will be refined using GPD tail.
% G1out  : Boolean indicating whether G1 should be removed from the null
%          distribution.
%
% Output:
% P      : P-values.
% apar   : Scale parameter of the GPD.
% kpar   : Shape parameter of the GPD.
% upar   : Location parameter of the GPD.
%
% For a complete description, see:
% * Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM.
%   Faster permutation inference in brain imaging.
%   Neuroimage. 2016 Jun 7;141:502-516.
%   http://dx.doi.org/10.1016/j.neuroimage.2016.05.068
% 
% This function is based on the following papers:
% * Knijnenburg TA, Wessels LFA, Reinders MJT, Shmulevich I. Fewer
%   permutations, more accurate P-values. Bioinformatics.
%   2009;25(12):i161-8.
% * Hosking JRM, Wallis JR. Parameter and quantile estimation for
%   the generalized Pareto distribution. Technometrics.
%   1987;29:339-349.
% * Grimshaw SD. Computing Maximum Likelihood Estimates for the
%   Generalized Pareto Distribution. Technometrics. 1993;35:185-191.
% * Choulakian V, Stephens MA. Goodness-of-Fit Tests for the
%   Generalized Pareto Distribution. Technometrics. 2001;43(4):478-484.
% 
% Also based on the tool released by Theo Knijnenburg, available at:
% https://sites.google.com/site/shmulevichgroup/people/theo-knijnenburg
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

% Compute the usual permutation p-values.
if G1out,
    Gdist = Gdist(2:end,:);
end
P    = palm_datapval(G,Gdist,rev);
Pidx = P < Pthr; % don't replace this "<" for "<=".

% If some of these are small (as specified by the user), these
% will be approximated via the GPD tail.
if any(Pidx),
    
    % Number of permutations & distribution CDF
    nP   = size(Gdist,1);
    if rev,
        [~,Gdist,Gcdf] = palm_competitive(Gdist,'descend',true);
    else
        [~,Gdist,Gcdf] = palm_competitive(Gdist,'ascend',true);
    end
    Gcdf = Gcdf/nP;
    
    % Keep adjusting until the fit is good. Change the step to 10 to get
    % the same result as Knijnenburg et al.
    Q  = (751:10:999)/1000;
    nQ = numel(Q);
    q  = 1;
    Ptail = NaN;
    while any(isnan(Ptail)) && q < nQ-1,

        % Get the tail
        qidx  = Gcdf >= Q(q);
        Gtail = Gdist(qidx);
        qi    = find(qidx,1);
        if qi == 1,
            upar = Gdist(qi) - mean(Gdist(qi:qi+1));
        else
            upar = mean(Gdist(qi-1:qi));
        end
        if rev,
            ytail = upar - Gtail;
            y     = upar - G((G < upar) & Pidx);
        else
            ytail = Gtail - upar;
            y     = G((G > upar) & Pidx) - upar;
        end
        
        % Estimate the distribution parameters. See §3.2 of Hosking &
        % Wallis (1987). Compared to the usual GPD parameterisation, 
        % here k = shape (xi), and a = scale.
        x    = mean(ytail);
        s2   = var(ytail);
        apar = x*(x^2/s2 + 1)/2;
        kpar =   (x^2/s2 - 1)/2;
        
        % Check if the fitness is good
        A2pval = andersondarling(gpdpvals(ytail,apar,kpar),kpar);
            
        % If yes, keep. If not, try again with the next quantile.
        if A2pval > .05;
            cte = numel(Gtail)/nP;
            Ptail = cte*gpdpvals(y,apar,kpar);
        else
            q = q + 1;
        end
    end
    
    % Replace the permutation p-value for the approximated
    % p-value
    if ~ isnan(Ptail),
        if rev,
            P((G < upar) & Pidx) = Ptail;
        else
            P((G > upar) & Pidx) = Ptail;
        end
    end
end

% ==============================================================
function p = gpdpvals(x,a,k)
% Compute the p-values for a GPD with parameters a (scale)
% and k (shape).
if abs(k) < eps;
    p = exp(-x/a);
else
    p = (1 - k*x/a).^(1/k);
end
if k > 0;
    p(x > a/k) = 0;
end

% ==============================================================
function A2pval = andersondarling(z,k)
% Compute the Anderson-Darling statistic and return an
% approximated p-value based on the tables provided in:
% * Choulakian V, Stephens M A. Goodness-of-Fit Tests
%   for the Generalized Pareto Distribution. Technometrics.
%   2001;43(4):478-484.

% This is Table 2 of the paper (for Case 3, in which 
% a and k are unknown, bold values only)
ktable = [0.9 0.5 0.2 0.1 0 -0.1 -0.2 -0.3 -0.4 -0.5]';
ptable = [0.5 0.25 0.1 0.05 0.025 0.01 0.005 0.001];
A2table = [ ...
    0.3390 0.4710 0.6410 0.7710 0.9050 1.0860 1.2260 1.5590
    0.3560 0.4990 0.6850 0.8300 0.9780 1.1800 1.3360 1.7070
    0.3760 0.5340 0.7410 0.9030 1.0690 1.2960 1.4710 1.8930
    0.3860 0.5500 0.7660 0.9350 1.1100 1.3480 1.5320 1.9660
    0.3970 0.5690 0.7960 0.9740 1.1580 1.4090 1.6030 2.0640
    0.4100 0.5910 0.8310 1.0200 1.2150 1.4810 1.6870 2.1760
    0.4260 0.6170 0.8730 1.0740 1.2830 1.5670 1.7880 2.3140
    0.4450 0.6490 0.9240 1.1400 1.3650 1.6720 1.9090 2.4750
    0.4680 0.6880 0.9850 1.2210 1.4650 1.7990 2.0580 2.6740
    0.4960 0.7350 1.0610 1.3210 1.5900 1.9580 2.2430 2.9220];

% The p-values are already sorted
k  = max(0.5,k);
z  = flipud(z)';
n  = numel(z);
j  = 1:n;

% Anderson-Darling statistic and p-value:
A2 = -n -(1/n)*((2*j-1)*(log(z) + log(1-z(n+1-j)))');
i1 = interp1(ktable,A2table,k,'linear','extrap');
i2 = interp1(i1,ptable,A2,'linear','extrap');
A2pval = max(min(i2,1),0);
