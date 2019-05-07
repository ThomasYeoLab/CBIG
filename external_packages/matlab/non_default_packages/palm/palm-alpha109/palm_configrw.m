function varargout = palm_configrw(varargin)
% Read and write PALM configuration files.
% 
% Usage:
% cfg = palm_configrw(fname)
% palm_configrw(cfg,fname)
% 
% cfg   : Configurations (cell array).
% fname : Text-file with the configurations.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Jan/2013
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

if nargin == 1,
    
    % Read files
    plmfile = varargin{1};
    fid = fopen(plmfile,'r');
    cfg = textscan(fid,'%s','CommentStyle','#');
    fclose(fid);
    if nargout == 1,
        varargout = cfg;
    end
    
elseif nargin == 2,
    
    % Version & environment
    ver = fliplr(strtok(fliplr(palm_help('version')),' '));
    ver = strtok(ver,')');
    if palm_isoctave,
        envrun = 'Octave';
    else
        envrun = 'MATLAB';
    end
    
    % Write files
    cfg = varargin{1};
    plmfile = varargin{2};
    fid = fopen(plmfile,'w');
    fprintf(fid,'# Configuration file for PALM.\n');
    fprintf(fid,'# Version %s, running in %s %s.\n',ver,envrun,version);
    fprintf(fid,'# %s\n',datestr(now));
    fprintf('Running PALM %s using %s %s with the following options:',ver,envrun,version);
    for c = 1:numel(cfg),
        s2d = str2double(cfg{c});
        if strcmp(cfg{c}(1),'-') && (isnan(s2d) || ~isreal(s2d)),
            fprintf(    '\n%s',cfg{c});
            fprintf(fid,'\n%s',cfg{c});
        else
            if ischar(cfg{c}),
                fprintf(    ' %s',cfg{c});
                fprintf(fid,' %s',cfg{c});
            else
                fprintf(    ' %s',num2str(cfg{c}));
                fprintf(fid,' %s',num2str(cfg{c}));
            end
        end
    end
    fprintf(    '\n');
    fprintf(fid,'\n');
    fclose(fid);
else
    error('Incorrect number of input arguments.');
end
