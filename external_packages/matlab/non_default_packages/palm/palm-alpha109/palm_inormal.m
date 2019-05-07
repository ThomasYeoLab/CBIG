function Z = palm_inormal(varargin)
% Applies a rank-based inverse normal transformation.
% 
% Usage: Z = inormal(X)
%            inormal(X,c)
%            inormal(X,method)
%            inormal(X,...,quanti)
%
% Inputs:
% X      : Original data. Can be a vector or an array.
% c      : Constant to be used in the transformation.
%          Default c=3/8 (Blom).
% method : Method to choose c. Accepted values are:
%              'Blom'   (c=3/8),
%              'Tukey'  (c=1/3),
%              'Bliss', (c=1/2)  and
%              'Waerden' or 'SOLAR' (c=0).
% quanti : All data guaranteed to be quantitative and
%          without NaN?
%          This can be a true/false. If true, the function
%          runs much faster if X is an array.
%          Default is false.
%
% Outputs:
% Z      : Transformed data.
% 
% References:
% * Van der Waerden BL. Order tests for the two-sample
%   problem and their power. Proc Koninklijke Nederlandse
%   Akademie van Wetenschappen. Ser A. 1952; 55:453-458
% * Blom G. Statistical estimates and transformed
%   beta-variables. Wiley, New York, 1958.
% * Tukey JW. The future of data analysis.
%   Ann Math Stat. 1962; 33:1-67.
% * Bliss CI. Statistics in biology. McGraw-Hill,
%   New York, 1967.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Aug/2011 (first version)
% Jun/2014 (this version)
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

% Accept inputs & defaults
c0 = 3/8;  % Default (Blom, 1958)
quanti = false;
if nargin == 1,
    c = c0;
elseif nargin > 1 && ischar(varargin{2}),
    switch lower(varargin{2}),
        case 'blom'
            c = 3/8;
        case 'tukey'
            c = 1/3;
        case 'bliss'
            c = 1/2;
        case 'waerden'
            c = 0;
        case 'solar'
            c = 0; % SOLAR is the same as Van der Waerden
        otherwise
            error('Method %s unknown. Use either ''Blom'', ''Tukey'', ''Bliss'', ''Waerden'' or ''SOLAR''.',varargin{2});
    end
elseif nargin > 1 && isscalar(varargin{2}),
    c = varargin{2};  % For a user-specified value for c
end
if nargin == 3,
    quanti = varargin{3};
    if isempty(varargin{2}),
        c = c0;
    end
end
X = varargin{1};

% If the trait is quantitative, avoid the loop
if quanti,
    % Get the rank for each value
    [~,iX] = sort(X);
    [~,ri] = sort(iX);
    
    % Do the actual transformation
    N = size(X,1);
    p = ((ri-c)/(N-2*c+1));
    Z = sqrt(2)*erfinv(2*p-1);
    
else
    Z = nan(size(X));
    for x = 1:size(X,2),
        
        % Remove NaNs
        XX = X(:,x);
        ynan = ~isnan(XX);
        XX = XX(ynan);
        
        % Get the rank for each value
        [~,iX] = sort(XX);
        [~,ri] = sort(iX);
        
        % Do the actual transformation
        N = size(XX,1);
        p = ((ri-c)/(N-2*c+1));
        Y = sqrt(2)*erfinv(2*p-1);
        
        % Check for repeated values
        [U,~,IC] = unique(XX);
        if numel(U) < N,
            sIC = sort(IC);
            dIC = diff(vertcat(sIC,1));
            U = unique(sIC(~dIC));
            for u = 1:numel(U),
                Y(IC == U(u)) = mean(Y(IC == U(u)));
            end
        end
        
        % Put the NaNs back
        Z(ynan,x) = Y;
    end
end
