function palm(varargin)
% Type 'palm' without arguments for help.
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
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

% This is redundant in when running as a function as all files should
% be together, but it helps when running from the shell
addpath(fileparts(mfilename('fullpath')));

% If Octave
if palm_isoctave,
    
    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);
    
    % If running as a script, take the input arguments
    [~,cmdname,~] = fileparts(program_invocation_name());
    if ~ strcmpi(cmdname(1:6),'octave'),
        varargin = argv();
    end
    
    % Be sure to print to the screen immediately
    page_screen_output(0);
    page_output_immediately(1);
    
    % Disable some warnings
    warning('off','Octave:precedence-change');
    warning('off','Octave:possible-matlab-short-circuit-operator');
    warning('off','Octave:function-name-clash');
    
else
    % This line marks the place up to nothing will be printed. It's long as
    % this because if it fails, at least it's not ugly and looks
    % purposeful from the outside :-)
    fprintf('.......................................................................\n');
end

% This is probably redundant but fixes a bug in an old Matlab version
nargin = numel(varargin);

% Print usage if no inputs are given
if nargin == 0 || strcmp(varargin{1},'-q'),
    palm_help('basic');
    return;
elseif nargin == 1,
    if any(strcmpi(varargin{1},{'-help','-?','-basic'})),
        palm_help('basic');
        return;
    elseif strcmpi(varargin{1},'-advanced'),
        palm_help('advanced');
        return;
    elseif strcmpi(varargin{1},'-checkprogs'),
        palm_checkprogs;
        return;
    end
end
palm_help('logo');

% Now run what matters
palm_core(varargin{:});
