function P = chi2Tests(X, method)
% Chi-square tests of homogeneity and independence
%   Computes the P-value for IxJ-table row/col independence.
%   Program by Steinar Thorvaldsen, steinar.thorvaldsen@uit.no, Dec. 2004. 
%   Ref.: DeltaProt toolbox at http://services.cbu.uib.no/software/deltaprot/
%   Last changes 22. Dec 2010.
%   Requires Matlab 7.1 or newer, and Matlab Statistics toolbox.
%   Input:
%        X:      data matrix (IxJ-table) of the observed frequency cells
%        method: 'RC': Read-Cressie power divergence statistics (default), lambda= 2/3 
%                'Pe': Standard Pearson chi2-distance, lambda= 1
%                'LL': Log Likelihood ratio distance, lambda= 0   
%   Output:
%        P-value
%
%   Use: P = chi2Tests(Observed,'RC') 

%	The P-value is computed through approximation with chi-2 distribution
%	under the null hypothesis for all methods.
%   The 'RC'-method is sligtly better than the 'Pe'-method in small tables 
%   with unbalanced column margins

%   Please, use the following reference:
%   Thorvaldsen, S. , Flå, T. and Willassen, N.P. (2010) DeltaProt: a software toolbox 
%       for comparative genomics. BMC Bioinformatics 2010, Vol 11:573.
%       See http://www.biomedcentral.com/1471-2105/11/573

%   Other reference:
%   Rudas, T. (1986): A Monte Carlo Comparision of Small Sample Behaviour
%       of The Pearson, the Likelihood Ratio and the Cressie-Read Statistics.
%       J.Statist. Comput. Simul, vol 24, pp 107-120.
%   Read, TRC and Cressie, NAC (1988): Goodness of Fit Statistics for
%       Discrete Multivariate Data. Springer Verlag.
%   Ewens, WJ and Grant, GR (2001): Statistical Methods in Bioinformatics.
%       Springer Verlag.

if nargin < 2
    method = 'RC'; %default value
else
    method=upper(method);
    switch method %validate the method parameter:
    case {'RC' 'PE' 'LL'}
        % these are ok
    otherwise
        error('Chi2Test:UnknownMethod', ...
          'The ''method'' parameter value must be ''RC'', ''Pe'', or ''LL''.');
    end %switch
end %if

% lambda - power for statistic:
if method == 'RC'
    lambda=2/3;
elseif method == 'PE'
    lambda=1;
elseif method == 'LL'
    lambda=0;
end

if any(any(X < 0))
   error('chi2Test expects counts that are nonnegative values');
   return;
end

[I J] = size(X);
if I < 2 | J < 2,
    error('Matrix of observation must at least be of size 2x2');
    return;
end

rs = sum(X')';
cs = sum(X);
Ns = sum(rs);
if Ns <= 0 % no counts found in table
    P=NaN;
    return;
end
Eo = (rs*cs)./Ns;  %matrix of null means

% Make sure expected values are not too small
Emin=1.5;
if any(any(Eo < Emin))
    %disp ('Note: Some expected values in chi2Test are small (<1.5)...');
end

if any(any(Eo <= 0))
   t=NaN;
elseif lambda == 0 % logL
   i = find(X>0);
   t = 2 * sum(X(i).*(log(X(i)./Eo(i)))); % the test statistic
elseif lambda == 1 % Pearson chi2
   t = sum(sum((Eo-X).^2 ./ Eo)); % Pearson test statistic
else % Read-Cressie test statistic
   t = (2/(lambda*(lambda+1)))*sum(sum(X.*((X./Eo).^lambda -1))); 
end

df=prod(size(X)-[1,1]); % degree of freedom
P = 1-chi2cdf(t,df);
