function g = palm_randg(a,siz)
% Generate random Gamma variables.
%
% g = palm_randg(a,siz)
%
% a   : Parameter A (= k) of the Gamma distribution. The shape
%       (or rate) parameter is fixed at 1.
% siz : A vector indicating the size of the array to be created.
% g   : Gamma-distributed variables.
%
% Reference:
% * Marsaglia G, Tsang WW. A simple method for generating gamma
%   variables. ACM Trans Math Softw. 2000;26(3):363-372.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% May/2015
% http://brainder.org

if exist('randg','builtin') || exist('randg','file'),
    
    % randg is a build-in in Octave, and a MEX in Matlab in the 
    % Stats Toolbox. If it's available, let's use it.  It's
    % about 4x faster in Octave, and 6x times faster in Matlab.
    g = randg(a,siz);
    
else
    
    % Constants for below
    if a < 1,
        d = 1+a-1/3;
    else
        d = a-1/3;
    end
    c = 1/sqrt(9*d);
    
    % Fill g as the random vars satisfy the criteria below
    g  = nan(siz);
    ig = isnan(g(:));
    while any(ig),
        
        gw = nan(sum(ig),1);
        
        % Make the Gaussian variable
        x  = nan(size(gw));
        v  = -ones(size(gw));
        iv = true(size(gw));
        while any(iv),
            x(iv) = randn(sum(iv),1);
            v(iv) = (1+c.*x(iv)).^3;
            iv(iv) = v(iv) <= 0;
        end
        
        % Make the uniform variable
        U = rand(size(gw));
        
        % Apply both criteria
        iU = U < 1-0.0331*x.^4;
        if any(iU),
            gw(iU) = d.*v(iU);
        end
        iU = ~iU;
        if any(iU),
            iU(iU) = log(U(iU)) < 0.5*x(iU).^2+d*(1-v(iU)+log(v(iU)));
            gw(iU) = d.*v(iU);
        end
        
        % Go for the next round if needed
        g(ig) = gw;
        ig    = isnan(g(:));
    end
    
    % See Marsaglia & Tsang (2000), page 371 (the note below the Table):
    % Gamma(a) = Gamma(a+1)*U^(1/a)
    if a < 1,
        g = g.*rand(size(g)).^(1/a);
    end
    
end
