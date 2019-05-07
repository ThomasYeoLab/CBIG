function P = palm_quicksave(X,flg,opts,plm,y,m,c,filename)
% An intermediate function to be used many times in palm.m,
% which chooses an adequate mask, then places back the data
% into the points defined by this mask, and save.
%
% Usage:
% P = palm_quicksave(X,flg,opts,plm,y,c,filename)
%
% Inputs:
% X        : Data to be saved.
% flg      : A flag that can be:
%            0: meaning that X is to be just saved.
%            1: meaning that X is a P-value, so it may
%               be converted to 1-P or to -log(P) if
%               this was specified by the user.
%            2: meaning that X is a statistic (G or Z)
%               that needs to be converted to a P-value,
%               and then, perhaps, to 1-p or -log(P).
% opts,plm : Structs with options and general data.
% y        : Index for the current data in plm.Yset.
% m        : Index for the current design in plm.Mset.
% c        : Index for the current contrast in plm.Cset.
% filename : File to be created.
%
% Outputs:
% P        : True P-value (not 1-p or -log(p)).
%
% _____________________________________
% Anderson M. Winkler
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

% Save only if X isn't empty
if ~ isempty(X),
        
    % Just some verbosity, if asked
    if opts.showprogress,
        fprintf('\t Saving file: %s\n',filename);
    end
    
    % Modify X accodring to the flag
    if flg == 1 
        
        % Deal with complex values in the gamma approximation
        if opts.accel.gamma,
            if ~ isreal(X),
                iidx  = imag(X(:)) ~= 0;
                nimag = sum(iidx);
                warning(sprintf([...
                    'Complex values found with the gamma approximation in %d out of %d tests. \n' ...
                    '         These will be saved as non-significant (p-value = 1).\n' ...
                    '         To solve this issue, increase the number of permutations.\n'...
                    '         Affected file: %s'],...
                    nimag,numel(iidx),filename));
                X(iidx) = 1;
                X(iidx) = real(X(iidx));
            end
        end
        
        % Convert P to 1-P
        if opts.savecdf,
            X = 1 - X;
        end
        
    elseif flg == 2,
        
        % Convert a statistic to a P-value
        X0 = X;
        if opts.savecdf,
            
            % CDF (1-P)
            X = palm_gcdf(X0,plm.rC{m}(c),plm.df2{y}{m}{c});
            
            % Even saving the CDF, the true p-vals may be needed
            if nargout > 0,
                P = palm_gpval(X0,plm.rC{m}(c),plm.df2{y}{m}{c});
            end
            
            % Double the p-vals for parametric 2-tailed
            if opts.twotail && plm.rC{m}(c) <= 1,
                upp     = X <= .5;
                X( upp) = 2*X;
                X(~upp) = 2*palm_gpval(X0,plm.rC{m}(c),plm.df2{y}{m}{c});
            end
        else
            % Just P
            X = palm_gpval(X0,plm.rC{m}(c),plm.df2{y}{m}{c});
            if nargout > 0,
                P = X;
            end
            
            % Double the p-vals for parametric 2-tailed
            if opts.twotail && plm.rC{m}(c) <= 1,
                upp     = X <= .5;
                X( upp) = 2*X;
                X(~upp) = 2*palm_gcdf(X0,plm.rC{m}(c),plm.df2{y}{m}{c});
            end
        end
    end
    
    % Convert to logarithm
    if opts.savelogp && any(flg == [1 2]),
        X = -log10(X);
    end
    
    % Prepare struct to save
    if opts.inputmv,
        S.filename = filename;
        S.readwith = 'load';
        S.data     = X;
    else
        % Choose an appropriate mask struct.
        if isempty(y), y = 1;  end
        if opts.npcmod || opts.MV || opts.CCA || opts.PLS || opts.forcemaskinter,
            S = plm.maskinter;
        else
            if plm.nmasks == 1,
                S = plm.masks{1};
            else
                S = plm.masks{y};
            end
        end
        
        % Inject the data.
        mask         = S.data;
        S.data       = double(S.data);
        S.data(mask) = X;
        if  strcmpi(lower(S.readwith),'dpxread') && (...
                strcmpi(plm.Ykindstr{y},'_dpv') || ...
                strcmpi(plm.Ykindstr{y},'_dpf') || ...
                strcmpi(plm.Ykindstr{y},'_dpx') ),
            S.filename = sprintf('%s.%s',filename,plm.Ykindstr{y}(2:end));
        else
            S.filename = filename;
        end
    end
    
    % Transpose back
    if opts.transposedata && max(size(S.data)) == numel(S.data),
        S.data = S.data';
    end
    
    % Save
    palm_miscwrite(S,true);
end