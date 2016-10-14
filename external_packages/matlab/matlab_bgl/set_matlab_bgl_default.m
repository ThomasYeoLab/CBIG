function old_default = set_matlab_bgl_default(varargin)
% SET_MATLAB_BGL_DEFAULT Sets a default option for the Matlab BGL interface
%
% old_default = set_matlab_bgl_default(options) or
% old_default = set_matlab_bgl_default(...) for key-value pair version
% options.istrans: the input matrices are already transposed [{0} | 1]
% options.nocheck: skip the input checking [{0} | 1]
% options.full2sparse: convert full matrices to sparse [{0} | 1]
%
% to get the current set of default options, call
% options = set_matlab_bgl_default()
%
% These options can make the Matlab BGL interface more efficient by
% eliminating the copying operations that occur between Matlab's structures
% and the BGL structures.  However, they are more difficult to use and are
% disabled by default.
%
% Generally, they are best used when you want to perform a large series of
% computations.
%
% Example:
%   % tranpose the matrix initially...
%   At = A'
%   old_options = set_matlab_bgl_default(struct('istrans',1));
%   % perform a bunch of graph work with At...
%   d1 = dfs(At,1); d2 = dfs(At,2); ...
%   % restore the old options 
%   set_matlab_bgl_default(old_options);

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%%

persistent default_options;
if ~isa(default_options,'struct')
    % initial default options
    default_options = struct('istrans', 0, 'nocheck', 0, 'full2sparse', 0);
end

if nargin == 0
    old_default = default_options;
else
    old_default = default_options;
    default_options = merge_options(default_options,varargin{:});
end
   
