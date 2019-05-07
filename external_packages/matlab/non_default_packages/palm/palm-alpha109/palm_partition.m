function [X,Z,eCm,eCx] = palm_partition(M,C,meth,Y)
% Partition a design matrix into regressors of interest and
% nuisance according to a given contrast.
% 
% Usage
% [X,Z] = palm_partition(M,C,meth,Y)
% 
% Inputs:
% M    : Design matrix, to be partitioned.
% C    : Contrast that will define the partitioning.
% meth : Method for the partitioning. It can be:
%        - 'Guttman'
%        - 'Beckmann'
%        - 'Winkler' (for testing only)
%        - 'Ridgway'
%        - 'none' (does nothing, X=M, Z=[])
% Y    : (Optional) For the 'Winkler' method only.
% 
% Outputs:
% X    : Matrix with regressors of interest.
% Z    : Matrix with regressors of no interest.
% eCm  : Effective contrast, equivalent to the original,
%        for the partitioned model [X Z], and considering
%        all regressors.
% eCx  : Same as above, but considering only X.
%
% References:
% * Guttman I. Linear Models: An Introduction. Wiley,
%   New York, 1982.
% * Smith SM, Jenkinson M, Beckmann C, Miller K,
%   Woolrich M. Meaningful design and contrast estimability
%   in FMRI. Neuroimage 2007;34(1):127-36.
% * Ridgway GR. Statistical analysis for longitudinal MR
%   imaging of dementia. PhD thesis. 2009.
% * Winkler AM, Ridgway GR, Webster MG, Smith SM, Nichols TE.
%   Permutation inference for the general linear model.
%   Neuroimage. 2014 May 15;92:381-97.
% _____________________________________
% A. Winkler, G. Ridgway & T. Nichols
% FMRIB / University of Oxford
% Mar/2012 (1st version)
% Aug/2013 (major revision)
% Dec/2015 (this version)
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

switch lower(meth),
    case 'guttman', % works for evperdat (3D)
        idx   = any(C~=0,2);
        X     = M(:,idx,:);
        Z     = M(:,~idx,:);
        eCm   = vertcat(C(idx,:),C(~idx,:));
        
    case 'beckmann', % works for evperdat (3D)
        Cu    = null(C');
        for t = 1:size(M,3),
            D        = pinv(M(:,:,t)'*M(:,:,t));
            CDCi     = pinv(C'*D*C);
            Pc       = C*CDCi*C'*D;
            Cv       = Cu - Pc*Cu;
            F3       = pinv(Cv'*D*Cv);
            if t == 1,
                X = zeros(size(M,1),size(CDCi,2),size(M,3));
                Z = zeros(size(M,1),size(F3,2),size(M,3));
            end
            X(:,:,t) = M(:,:,t)*D*C*CDCi;
            Z(:,:,t) = M(:,:,t)*D*Cv*F3;
        end
        eCm = vertcat(eye(size(X,2)),...
            zeros(size(Z,2),size(X,2)));
        
    case 'winkler', % this is just for testing, and doesn't work with evperdat
        D     = pinv(M'*M);
        X     = M*D*C*pinv(C'*D*C);
        Z     = (M*D*M'-X*pinv(X))*Y;
        eCm   = vertcat(eye(size(X,2)),...
            zeros(size(Z,2),size(X,2)));
        
    case 'ridgway', % works for evperdat (3D)
        rC    = rank(C);
        rM    = rank(M(:,:,ceil(size(M,3)/2)));
        rZ    = rM - rC;
        pinvC = pinv(C');
        C0    = eye(size(M,2)) - C*pinv(C);
        for t = 1:size(M,3),
            if t == 1,
                X = zeros(size(M,1),size(pinvC,2),size(M,3));
                Z = zeros(size(M,1),rZ,size(M,3));
            end
            tmpX = M(:,:,t)*pinvC;
            [tmpZ,~,~]  = svd(M(:,:,t)*C0);
            Z(:,:,t) = tmpZ(:,1:rZ);
            X(:,:,t) = tmpX-Z(:,:,t)*pinv(Z(:,:,t))*tmpX;
        end
        eCm = vertcat(eye(size(X,2)),...
            zeros(size(Z,2),size(X,2)));
        
    case 'none', % works for evperdat (3D)
        X     = M;
        Z     = [];
        eCm   = C;
        
    otherwise
        error('''%s'' - Unknown partitioning scheme',meth);
end
eCx = eCm(1:size(X,2),:);
